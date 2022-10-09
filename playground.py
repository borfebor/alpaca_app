#needed packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from random import randint
from math import pi
import seaborn as sns
from src.alpaca_new import alpaca
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
    df = alpaca.eats(uploaded_file)
     
    st.sidebar.header('Data preprocessing')
    cleaning = st.sidebar.checkbox('Data cleaning', value=True)
    formating = st.sidebar.checkbox('Data formating', value=True)
    
    lfq_method = st.sidebar.selectbox('Label-Free Quantification method', ['iBAQ', 'LFQ'], 0)
    
    sub = st.sidebar.empty()
    
    subpro_approach = sub.checkbox('Subproteomic analysis', False)
    if subpro_approach == True:
          subproteome = sub.text_input('Analysed subproteome')
    else:
          subproteome = None
    
    df, condition = alpaca.spits(df, lfq_method=lfq_method, cleaning=cleaning, formating=formating, subproteome=subproteome)
    st.title('Your data')
    st.write(df)
    st.write(f'Your data contains {len(condition)} experimental conditions ({condition})')
    
    experimenter = st.expander('Experimental set-up', expanded=True)
    with experimenter:
          st.subheader('Quantification standards parameters')
          replicate_list = list(df.Replicate.unique())
          replicate_list.append('All')
          replicate = st.selectbox('Were the quantification standards added to just a replicate per set?', replicate_list, 0)
          concentration = st.number_input('Quantification Standards concentration', 0.0, 21.2, 0.5)
          volume = 10.6 / concentration
          
          amount = st.number_input('How much standard volume - in µl - have you spiked in your samples?', 0.0, volume, 6.0)
          
          st.header("Enrichment preparation")
          enrichments = st.number_input('Different sample preparations', 1, 10, 1)
          
          enrichment_t = [f'Preparation_{en+1}' for en in range(enrichments)]
         
          c1, c3, c4, c5 = st.columns([3,1,1,1])
          adder = 0
          enrichment_type_dict = dict()
          if subproteome != None:
               subproteome_dict = dict()
          else:
               subproteome_dict = None
          prep = list()
          for c in range(enrichments):
               enrich = c1.multiselect(enrichment_t[c], condition, condition[adder])
               dil = c3.number_input('Dilution', 1, 1000, 10 + adder)
               vol = c4.number_input('Spiked standard volume (µl)', 0.0, 1000.0, 8.5 + adder)
               sampl_v = c5.number_input('Sample volume (ml)', 0.0, 100.0, 45.0 + adder)
               prepa = [dil, vol, sampl_v]
               prep.append(prepa)
               enrichment_type_dict[enrichment_t[c]] = enrich
               if subproteome != None:
                    my_subproteome = list(df.Subproteome.unique())
                    subprot = st.selectbox('Prepared subproteome', my_subproteome, 0)
                    subproteome_dict[enrichment_t[c]] = subprot
               adder += 1
          
          sample_preparation = alpaca.scientist(prep, enrichment_type_dict, subproteome_dict)
          st.write(sample_preparation)
          
     expander1 = st.expander('Quantification')
     expander2 = st.expander('Enrichment')
     expander3 = st.expander('Absolute Protein Quantification')
     with expander1:
          st.header('Quantification')
     #UPS parameters input
          ups2 = pd.read_csv('ups2_dynamicStd.csv', sep=',')
     #UPS calculations based on the spiked amount of standards
          if replicate == 'All':
               replicate = None
          ups2, data, coef, inter = alpaca.census(df, 'ups2_dynamicStd.csv', concentration=concentration, 
                                                 in_sample=amount, lfq_col=lfq_method, replicate=replicate, save=False)

          st.write(data)
    
