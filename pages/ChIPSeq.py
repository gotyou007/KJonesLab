from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from pandas import read_csv,DataFrame
from rpy2.robjects.packages import importr # Loading R library, equivalent to library() in R
import streamlit as st
from rpy2.robjects.conversion import localconverter
from .utils import get_data_from_excel, convert_df
import plotly.express as px
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
    @st.cache(suppress_st_warning=True, allow_output_mutation=True)
    def process_files(bedfile ='VCP.bed'):
        with localconverter(ro.default_converter):
            CSK = importr('ChIPseeker')
            base = importr('base')
            GF = importr('GenomicFeatures')
            txdb = GF.makeTxDbFromGFF('gencode.v36.annotation.gff3.gz')
            meth = importr('methods')
            peaks = CSK.readPeakFile(bedfile)
            annoPeakTab = base.as_data_frame(meth.slot(CSK.annotatePeak(peaks,TxDb = txdb),"anno"))
        with localconverter(ro.default_converter + pandas2ri.converter):
            annoPeakTab = ro.conversion.rpy2py(annoPeakTab)
        return annoPeakTab
    annoPeakTab = process_files('VCP.bed')
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
    anno = annoPeakTab['annotation'].str.replace(r"Intron\s\(.*\)","Intron")
    anno = anno.str.replace(r"Exon\s\(.*\)","Exon")
    
    c = Counter(anno)
    fig1 = px.pie(values = [float(v) for v in c.values()], names = [k for k in c])
    
    
    #plot distance to TSS
    annoPeakTab['distSummary'] = pd.cut(
        x=annoPeakTab["distanceToTSS"],
        bins=[-np.inf, -100000,-10000,-3000,-1000,0,1000, 3000, 10000,100000, np.inf],
        labels=["<100 kb", "-100:-10 kb", "-10:-3 kb", '-3:-1kb','-1:0 kb', '0:1 kb', '1:3 kb', '3:10 kb', '10:100 kb', ">100 kb"]
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
    '''uploaded_files = st.file_uploader("Choose input file", accept_multiple_files=True)
    for uploaded_file in uploaded_files:
        bytes_data = uploaded_file.read()
        st.write("filename:", uploaded_file.name)
        st.write(bytes_data)'''
    ####use bigwig to compute socre
    def scoreBedByBigwig(bed, bigwig):
        bw = pyBigWig.open(bigwig)
        df = pd.read_csv("ChIP/Peaks/V5.bed", sep="\t", header=None)
        df.iloc[:,4] = df.apply(lambda row: bw.stats(row[0], int(row[1]), int(row[2]))[0], axis = 1)
        
        return df
    #upload bigwig and bed
    files = st.file_uploader("Make sure you upload bed file first, then bigwig", accept_multiple_files=True)

    if len(files) == 0:
        st.error("No file were uploaded")

    for i in range(len(files)):
        bytes_data = files[i].getvalue()  # read the content of the file in binary
    
        with open(os.path.join("/tmp", files[i].name), "wb") as f:
            f.write(bytes_data)  # write this content elsewhere
    with st.spinner('Wait for it...'):
        if st.button('Run computeMatrix'):
            outputMatrix='tempMatrix.gz'
            cmd1 = ['computeMatrix', 'reference-point', '--referencePoint', 'center', '-b', '4000', '-a', '4000', 
                    '-R', os.path.join("/tmp", files[0].name), '-S', os.path.join("/tmp", files[1].name), '--sortRegions', 'descend', '--skipZeros', '-o', 'tempMatrix.gz',
                    '--outFileSortedRegions', 'sorted.bed']
            subprocess.run(cmd1)
            if os.path.isfile('tempMatrix.gz'):
                cmd1 = ['plotHeatmap', '-m', 'tempMatrix.gz', '-out', 'tornadoPlot.eps', '--sortRegions', 'keep', '--colorList', 'white,orange', 
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


    
    
    # ---- HIDE STREAMLIT STYLE ----
    hide_st_style = """
                <style>
                #MainMenu {visibility: hidden;}
                footer {visibility: hidden;}
                header {visibility: hidden;}
                </style>
                """
    st.markdown(hide_st_style, unsafe_allow_html=True)

