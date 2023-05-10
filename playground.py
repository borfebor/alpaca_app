#needed packages
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image

from src.alpaca_new import alpaca
from src.viz import Viz


st.set_page_config(
     page_title="Alpaca Proteomics",
     layout="wide",
     initial_sidebar_state="expanded"
)

image = Image.open('ALPACA_LOGO2.png')
icon = Image.open('ALPACA_ICON.png')
alpaca_enrichment = Image.open('enrichment_instructions.png')
instructions = Image.open('instructions.png')
st.sidebar.image(image)


top_bar = st.sidebar.empty()
uploaded_file = top_bar.file_uploader('Import your Protein Groups file:',
                                        type=['csv','txt','xlsx'],
                                         help='You can upload directly your ProteinGroup.txt from MaxQuant output '
                                              'and clean it using the button below, or pre-clean the data with Perseus. '
                                              'Our quantification approach uses iBAQ for Absolute Quantification.')

example_data = None

if uploaded_file == None:
    
    example_data = st.sidebar.selectbox('Example datasets', [None, 'Ferrero-Bordera et al. 2023', 'Maaß et al. 2014'])
     
if uploaded_file is not None:

    # Can be used wherever a "file-like" object is accepted:
    df = alpaca.eats(uploaded_file)
     
    st.sidebar.header('Data preprocessing')
    cleaning = st.sidebar.checkbox('Data cleaning', value=True)
    formatting = st.sidebar.checkbox('Data formatting', value=True)
    
    lfq_options = ['iBAQ', 'LFQ', 'Top3', 'Intensity', 'MS/MS count']
    
    lfq_menu = list(dict.fromkeys([item for item in lfq_options for col in df.select_dtypes(include=np.number).columns if item in col]))
    
    lfq_method = st.sidebar.selectbox('Label-Free Quantification method', lfq_menu, 0)
    
    normalization = st.sidebar.selectbox('Intensity normalization', [False, 'Relative', 'Median', 'Quantile'], 0)
      
    st.write(df)
    df, condition, lfq_method = alpaca.spits(df, lfq_method=lfq_method, cleaning=cleaning, normalization=normalization,
                                formatting=formatting, identifier=None)

