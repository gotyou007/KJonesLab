from optparse import Values
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from pandas import read_csv,DataFrame
from rpy2.robjects.packages import importr # Loading R library, equivalent to library() in R
import streamlit as st
from rpy2.robjects.conversion import localconverter
from .utils import get_data_from_excel, convert_df, convert_df_noindex
import plotly.express as px
import plotly.graph_objects as go
from collections import Counter
import pandas as pd
import numpy as np
import subprocess
import pyBigWig
import os
import tempfile
from pathlib import Path




def app():
    """
    ChIP peaks analysis 
    """
    
    ##input bed file for annotation
    uploaded_file = st.file_uploader("Choose a file for peaks annotation analysis")
    st.sidebar.header("Please pick species:")
    species = st.sidebar.selectbox('Species',('human', 'mouse'))
    st.sidebar.header("Please select parameters for your analysis:")
    peaks_type = st.sidebar.selectbox('peaks type',("Promoter", "Distal Intergenic"))
    referencePoint = st.sidebar.selectbox('Reference Point?',("TSS", "center", "TES"))
    sorted_output = st.sidebar.selectbox('sort bed?',("Yes", "No"))
    tornado_color = st.sidebar.selectbox('tornado plot color',('white,red', 'white,orange', 'white,green', 'white,purple', 'white,yellow'))
    @st.cache(suppress_st_warning=True, allow_output_mutation=True)
    def prepare_anno(species):
        GFF = 'gencode.vM10.annotation.gff3.gz' if species == 'mouse' else "gencode.v36.annotation.gff3.gz"
        with localconverter(ro.default_converter):
            CSK = importr('ChIPseeker')
            base = importr('base')
            GF = importr('GenomicFeatures')
            txdb = GF.makeTxDbFromGFF(GFF)
            meth = importr('methods')
            
        return CSK, base, GF, txdb, meth
    def process_files(bedfile, species):
        CSK, base, GF, txdb, meth = prepare_anno(species)
        with localconverter(ro.default_converter):
            peaks = CSK.readPeakFile(bedfile)
            annoPeakTab = base.as_data_frame(meth.slot(CSK.annotatePeak(peaks,TxDb = txdb),"anno"))
        with localconverter(ro.default_converter + pandas2ri.converter):
            annoPeakTab = ro.conversion.rpy2py(annoPeakTab)
        return annoPeakTab
    def extract_peaks_by_type(peaks: object, anno: str) -> object:
        df = peaks.loc[peaks['annotation'] == anno].iloc[:,range(3)]
        return df

    with st.spinner('Wait for it...'):
        if st.button('Run peaks annotation process'):
            if uploaded_file is not None:
                bytes_data = uploaded_file.getvalue()  # read the content of the file in binary
            
                with open(os.path.join("/tmp", uploaded_file.name), "wb") as f:
                    f.write(bytes_data)  # write
        
            annoPeakTab = process_files(os.path.join("/tmp", uploaded_file.name), species)
            st.dataframe(annoPeakTab[:5])
            csv = convert_df(annoPeakTab)
            
            st.download_button(
            "Press to Download peaks annoTab",
            csv,
            file_name="peakAnno.csv",
            mime='text/csv',
            key='download-csv')
    
            #visualize anno
            annoPeakTab['annotation'] = annoPeakTab['annotation'].astype('string')
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Intron\s\(.*\)","Intron")
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Exon\s\(.*\)","Exon")
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Promoter\s\(.*\)","Promoter")
            c = Counter(annoPeakTab['annotation'])

            #need a specific type of peaks (e.g. promoter, itergenic ...)
            s_peaks = extract_peaks_by_type(annoPeakTab, peaks_type)
            s_peaks = convert_df_noindex(s_peaks)
            st.download_button(
            "Press to Download your type-peaks annoTab",
            s_peaks,
            file_name="type_peakAnno.csv",
            mime='text/csv',
            key='download-csv')
            
            keys = ['Promoter', 'Distal Intergenic', 'Intron', 'Exon', "3' UTR", 'Downstream (<=300bp)',"5' UTR"]
            values = [float(c[k]) for k in keys]
            colors = px.colors.qualitative.Plotly[:7]
            
            fig1 = go.Figure(data=[go.Pie(labels=keys, values=values)])
            fig1.update_traces(marker=dict(colors=colors))
        

            #plot distance to TSS
            annoPeakTab['distSummary'] = pd.cut(
                x=annoPeakTab["distanceToTSS"],
                bins=[-np.inf, -100000,-10000,-3000,-1000,0,1000, 3000, 10000,100000, np.inf],
                labels=["<-100 kb", "-100:-10 kb", "-10:-3 kb", '-3:-1kb','-1:0 kb', '0:1 kb', '1:3 kb', '3:10 kb', '10:100 kb', ">100 kb"]
            )
    
            totalPeakNumber = annoPeakTab.shape[0]
            distSum = annoPeakTab.groupby(['distSummary']).size().reset_index()
            distSum['percentage'] = annoPeakTab.groupby(['distSummary']).size().apply(lambda x: 100 * x / float(totalPeakNumber)).values
            distSum.columns = ['distCategory', 'Count', 'Percentage']
            fig2 = px.bar(distSum, x='Count', color='distCategory',y=['s' for _ in range(distSum.shape[0])],text=distSum['Percentage'].apply(lambda x: '{0:1.2f}%'.format(x)),
                        orientation='h', height=400)
            fig2.update_layout(yaxis={'visible': False, 'showticklabels': False},
                            plot_bgcolor="rgba(0,0,0,0)",
                            xaxis = dict(
                            title = '',
                            tickmode = 'array',
                            tickvals = [0, sum(distSum.Count[:5]), distSum.Count.sum()],
                            ticktext = ['Upstream', 'TSS', 'Downstream'])
                            )
            st.dataframe(distSum)
            st.plotly_chart(fig1, use_container_width=True)
            st.plotly_chart(fig2, use_container_width=True)
################################################################
##deeptools::computing score on a bed file and make tornado plot
################################################################
    ####use bigwig to compute socre
    
    #upload bigwig and bed
    files = st.file_uploader("Make sure you upload bed file first, then bigwig (choose qvalue bigwig if you want to get score, choose fragment if you are plotting tornado", accept_multiple_files=True)

    if len(files) == 0:
        st.error("No file were uploaded")

    for i in range(len(files)):
        bytes_data = files[i].getvalue()  # read the content of the file in binary
    
        with open(os.path.join("/tmp", files[i].name), "wb") as f:
            f.write(bytes_data)  # write this content elsewhere

    with st.spinner('Wait for it...'):
        if st.button('Run computeMatrix'):
            outputMatrix='tempMatrix.gz'
            if sorted_output == 'Yes':
                cmd1 = ['computeMatrix', 'reference-point', '--referencePoint', referencePoint, '-b', '10000', '-a', '10000', 
                        '-R', os.path.join("/tmp", files[0].name), '-S', os.path.join("/tmp", files[1].name), '--sortRegions', 'descend', '--skipZeros', '-o', 'tempMatrix.gz',
                        '--outFileSortedRegions', 'sorted.bed']
            else:
                cmd1 = ['computeMatrix', 'reference-point', '--referencePoint', referencePoint, '-b', '10000', '-a', '10000', 
                        '-R', os.path.join("/tmp", files[0].name), '-S', os.path.join("/tmp", files[1].name), '--sortRegions', 'keep', '--skipZeros', '-o', 'tempMatrix.gz'
                        ]
            subprocess.run(cmd1)
            if os.path.isfile('tempMatrix.gz'):
                cmd1 = ['plotHeatmap', '-m', 'tempMatrix.gz', '-out', 'tornadoPlot.eps', '--sortRegions', 'keep', '--colorList', tornado_color, 
                        '--missingDataColor', '1', '--legendLocation', 'none', '--plotFileFormat', 'eps']

                subprocess.run(cmd1)
            st.success('Done!')
    if os.path.isfile('tornadoPlot.eps'):
        with open("tornadoPlot.eps", "rb") as file:
            btn = st.download_button(
                label="Download image",
                data=file,
                file_name="tornadoHeatmap.eps",
                mime="image/png"
            )
    def scoreBedByBigwig(bed, bigwig):
        bw = pyBigWig.open(bigwig)
        df = pd.read_csv(bed, sep="\t", header=None)
        scores = df.apply(lambda row: bw.stats(row[0], int(row[1]), int(row[2]))[0], axis = 1)
        
        return list(scores)
    
    #plot peaks in distance and qvalue
    with st.spinner('Wait for it...'):
        if st.button('Run scatter plotting for peaks'):
            scores = scoreBedByBigwig(os.path.join("/tmp", files[0].name), os.path.join("/tmp", files[1].name))
            annoPeakTab = process_files(os.path.join("/tmp", files[0].name), species)
            annoPeakTab['score'] = scores
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Intron\s\(.*\)","Intron")
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Exon\s\(.*\)","Exon")
            annoPeakTab['annotation'] = annoPeakTab['annotation'].str.replace(r"Promoter\s\(.*\)","Promoter")
            annoPeakTab['dist'] = annoPeakTab['distanceToTSS'].abs()+1

            fig3 = px.scatter(annoPeakTab, x="V5", y="dist", color='annotation',hover_data=['V4'])
            fig3.update_yaxes(type="log")
            fig3.update_traces(marker={'size': 3})
            
            fig4 = px.histogram(annoPeakTab, y="dist", color="annotation", barmode="overlay")
            fig4.update_layout(xaxis_range=[0,600])
            #fig.update_yaxes(type="log")
            
            
            distSum = annoPeakTab.groupby(['annotation']).size().reset_index(name='count')
            st.dataframe(distSum)
            st.plotly_chart(fig3, use_container_width=True)
            st.plotly_chart(fig4, use_container_width=True)



    
    
    # ---- HIDE STREAMLIT STYLE ----
    hide_st_style = """
                <style>
                #MainMenu {visibility: hidden;}
                footer {visibility: hidden;}
                header {visibility: hidden;}
                </style>
                """
    st.markdown(hide_st_style, unsafe_allow_html=True)

