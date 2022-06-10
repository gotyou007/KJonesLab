# @Email:  contact@pythonandvba.com
# @Website:  https://pythonandvba.com
# @YouTube:  https://youtube.com/c/CodingIsFun
# @Project:  Sales Dashboard w/ Streamlit



import pandas as pd  # pip install pandas openpyxl
import plotly.express as px  # pip install plotly-express
import streamlit as st  # pip install streamlit
import numpy as np

# emojis: https://www.webfx.com/tools/emoji-cheat-sheet/
st.set_page_config(page_title="DeSeq", page_icon=":bar_chart:", layout="wide")

# ---- READ EXCEL ----
@st.cache
def get_data_from_excel():
    df = pd.read_excel(
        io="pooled_DESeq.xlsx",
        engine="openpyxl",
        sheet_name="comb_vs_scr",
        usecols="A:K"
    )
    return df

df = get_data_from_excel()
# ---- SIDEBAR ----
st.sidebar.header("Please Filter Here:")
pvalue = st.sidebar.selectbox('P.adj',(0.01,0.05))


df_selection = df.query(
    "padj <= @pvalue"
)
df_selection['logpadj'] = -np.log10(df_selection['padj'])
# ---- MAINPAGE ----
st.title(":Volcano Plot")
st.markdown("##")
fig = px.scatter(df_selection, x="log2FoldChange", y="logpadj",
                     width=1200, height=800,
                     hover_name='gene_name')
                     
fig.update_layout(
    xaxis=dict(tickmode="linear"),
    plot_bgcolor="rgba(0,0,0,0)",
    yaxis=(dict(showgrid=False)),
    xaxis_range=[-4,4],
)

st.plotly_chart(fig, use_container_width=True)

st.markdown("""---""")



# ---- HIDE STREAMLIT STYLE ----
hide_st_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_st_style, unsafe_allow_html=True)
