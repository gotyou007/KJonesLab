import streamlit as st
import numpy as np
import pandas as pd
import seaborn as sns
import os
import plotly.express as px  # pip install plotly-express
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
from matplotlib.colors import ListedColormap
from scipy.stats import zscore
from sklearn import cluster

def app():
    def get_data_from_excel(sheet):
        df = pd.read_excel(
            io="pooled_DESeq.xlsx",
            engine="openpyxl",
            sheet_name=sheet
        )
        return df
    sheet = st.text_input("select sheet", 'comb_vs_scr')
    data = get_data_from_excel(sheet)
    st.sidebar.header("Please Filter Here:")
    pvalue = st.sidebar.selectbox('P.adj',(0.001,0.01,0.05))
    lfc = st.sidebar.selectbox('lfc',(0.5, 1, 2))
    label=("up", "not significant", "down")
    data.loc[(data["log2FoldChange"] >= lfc) & (data['padj'] < pvalue), 'label'] = label[0]  # upregulated
    data.loc[(data["log2FoldChange"] <= -lfc) & (data['padj'] < pvalue), 'label'] = label[2]  # downregulated
    data['label'].fillna(label[1], inplace=True)  # intermediate
    #export gene list
    dif_genes = data.loc[data['label'] != label[1]]

    #df_selection = data.query(
    #    "padj <= @pvalue"
    #)
    data['-logpadj'] = -np.log10(data['padj'])
    # ---- MAINPAGE ----
    st.title("Volcano Plot")
    st.markdown("####")
    
    color=("red", "grey", "green")
    cmap = dict(zip(label, color))
    fig = px.scatter(data, x="log2FoldChange", y="-logpadj",
                        width=1200, height=800,
                        hover_name='gene_name',
                        color='label', color_discrete_map=cmap)
                        
    fig.update_layout(
        xaxis=dict(tickmode="linear"),
        plot_bgcolor="rgba(0,0,0,0)",
        yaxis=(dict(showgrid=False)),
        xaxis_range=[-4,4],
    )
    @st.cache
    def convert_df(df):
        return df.to_csv().encode('utf-8')
    csv = convert_df(dif_genes)

    st.download_button(
    "Press to Download differential expressed genes",
    csv,
    "file.csv",
    "text/csv",
    key='download-csv'
    )

    st.plotly_chart(fig, use_container_width=True)

    st.markdown("""---""")

    #heatmap on kmean clustering
    #require df format: rows represent genes with gene name as index, and columns samples, values are rlog
    #sample_info:id, name, trt
    sample_info = get_data_from_excel('samples')
    
    sheet_2 = st.text_input("select sheet for heatmap plot", 'rlog')
    df_heat = get_data_from_excel(sheet_2)
    df_heat = df_heat.drop(columns=['id', 'biotype'])
    df_heat = df_heat.set_index('gene_name') #concert 'gene_name' to index
    df_heat = zscore(df_heat, axis=1) #z norm the df
    #convert column names from id to biological name
    df_heat.rename(columns = dict(zip(sample_info.id.to_list(), sample_info.name.to_list())), inplace=True)


    # Perform k-means clustering by using a pre-defined number of clusters
    k = st.slider('How many cluster you wish?', 3, 8, 5)
    def kmeans(df, k=5):
        '''
        take cleaned df and k number, perform k mean clustering on
        genes, with expression level in samples as features 
        '''
        kmeans = cluster.KMeans(n_clusters=k)
        kmeans.fit(df)
        # Extract labels and centroids
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        df['cluster_k'] = labels
        df_kmean = df.sort_values(by=['cluster_k'])
        return df_kmean
    
    df_heat = kmeans(df_heat, k)
    st.dataframe(sample_info[:5])

    #color map for kmean
    def color_mapping(obj, color_set):
        color = sns.color_palette("Paired", obj.unique().size)
        lut = dict(zip(obj.unique(), color))
        colors = obj.map(lut)
        return lut, colors
        

    kcluster = df_heat.pop('cluster_k') #we don't need cluster_k in plot, so we use pop function to fetch cluster information
    
    r_lut, row_colors = color_mapping(kcluster, "hls")
    #color map for samples
    sample_info = sample_info.iloc[:,1:3] #column 4 is factor size
    sample_info = sample_info.set_index('name') #use bilogical name as index
    trt = sample_info.trt
    #trt_color = sns.color_palette("Paired", trt.unique().size)
    #c_lut = dict(zip(trt.unique(), trt_color))
    #col_colors = trt.map(c_lut)
    c_lut, col_colors = color_mapping(trt, "Paired")

    #sns.clustermap(df_heat, row_cluster=False, cmap='BrBG', 
    #              row_colors=row_colors, col_colors=col_colors)
    #st.pyplot()
    st.set_option('deprecation.showPyplotGlobalUse', False)

    g = sns.clustermap(df_heat, row_cluster=False, cmap='bwr', 
                   row_colors=row_colors, col_colors=col_colors)
    for label in kcluster.unique():
        g.ax_col_dendrogram.bar(0, 0, color=r_lut[label],
                                label=label, linewidth=0)
    l1 = g.ax_col_dendrogram.legend(title = "cluster", ncol=1,
                                bbox_to_anchor=(0.1, 0.6), bbox_transform=gcf().transFigure)
    for label in trt.unique():
        g.ax_row_dendrogram.bar(0, 0, color=c_lut[label], 
                                label=label, linewidth=0)

    l2 = g.ax_row_dendrogram.legend(title = "trt", ncol=1, bbox_to_anchor=(0.1, 0.8), bbox_transform=gcf().transFigure)
    st.pyplot()



    # ---- HIDE STREAMLIT STYLE ----
    hide_st_style = """
                <style>
                #MainMenu {visibility: hidden;}
                footer {visibility: hidden;}
                header {visibility: hidden;}
                </style>
                """
    st.markdown(hide_st_style, unsafe_allow_html=True)