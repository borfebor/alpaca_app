#needed packages
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from random import randint
from math import pi
#import seaborn as sns
from src.alpaca import alpaca
from st_aggrid import AgGrid
import altair as alt
import time

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

uploaded_file = st.sidebar.file_uploader('Import your Protein Groups file:',
                                        type=['csv','txt','xlsx'],
                                         help='You can upload directly your ProteinGroup.txt from MaxQuant output '
                                              'and clean it using the button below, or pre-clean the data with Perseus. '
                                              'Our quantification approach uses iBAQ for Absolute Quantification.')

if uploaded_file is None:
    st.image(instructions)
    st.stop()

if uploaded_file is not None:

    # Can be used wherever a "file-like" object is accepted:
    df = alpaca.importer(uploaded_file)
    df, columns, default = alpaca.formater(df)

    n, r, condition = alpaca.experimenter(df)
    experimenter = st.expander('Experimental set-up', expanded=True)
    
    col = st.sidebar.multiselect('Choose the columns to import', columns, default)
    df = df[col]
     
    with experimenter:
        
         st.subheader('UPS2')
         volume = st.number_input('How much volume - in µl - have you used to resuspend your UPS?', 0.0, 50.0, 21.2)
         amount = st.number_input('How much standard volume - in µl - have you spiked in your samples?', 0.0, volume, 6.0)
        
         st.subheader('Cell count')
         c1, c2 = st.columns(2)
         default_count = 1060000000
         adder = 0
         condition_count = list()
         cell_count = list()
         for c in condition:
               cond = c1.multiselect('Counted condition', condition, c)
               cells = c2.number_input('Cell count (cells/ml)', 10000, 1000000000000, default_count + adder)
               cell_count.append(cells)
               condition_count = condition_count + cond
               adder += 10
          count_dict = dict(zip(condition_count, cell_count))
    
