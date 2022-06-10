
#needed packages

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from random import randint
from math import pi
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
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
    col = st.sidebar.multiselect('Choose the columns to import', columns, default)
    df = df[col]

cleaning = st.expander('Data Cleaning & Manipulation', expanded=True)
about = st.sidebar.expander('ALPACA', expanded=True)
with about:
    st.subheader('About Alpaca')
    st.markdown('Alpaca is coded with all my heart to ease your data analysis of '
                'Absolute Protein Quantification experiments. '
                'Thank you for pushing Systems Biology a bit further with your data.')

    st.markdown('[**Borja Ferrero-Bordera**](https://www.linkedin.com/in/borjaferrero/) 👨‍🔬🧬')
with cleaning:
    st.subheader('Data Cleaning & Manipulation')
    st.write('Here you can arrange your data. You can clean it out of Possible Contaminants, Reverses and Decoys, or '
             'transform it to a logaritmic scale for a proper linearity of your quantification.')
    data_cleaning = st.checkbox('Remove Possible Contaminants, Reverses and Decoys', False)
    log_transform = st.checkbox('Log2 transform iBAQ intensities', False)
    if st.button('Arrange my data for me 👍'):
        # Lines for data cleaning based in columns
        my_bar = st.empty()
        my_bar.progress(0)
        for percent_complete in range(100):
            time.sleep(0.01)
            my_bar.progress(percent_complete + 1)
        if data_cleaning == True & log_transform == True:
            df = alpaca.data_cleaner(df)
            df = alpaca.log_transform(df)
        elif data_cleaning == True:
            df = alpaca.data_cleaner(df)
        elif log_transform == True:
            df = alpaca.log_transform(df)

st.subheader('Data Overview')
st.write('Data Dimension: ' + str(df.shape[0]) + ' rows and ' + str(df.shape[1]) + ' columns.')
placeholder = st.empty()
placeholder.write(df)
n, r, condition = alpaca.experimenter(df)
rep = alpaca.replicator(r)
st.write('Your data contains ' + str(n) + ' experimental conditions with ' + str(r) + ' replicates each: \n' + str(condition))

expander1 = st.expander('Quantification')
expander2 = st.expander('Enrichment')
expander3 = st.expander('Absolute Protein Quantification')
expander4 = st.expander('Data visualization')

with expander1:
    st.header('Quantification')
#st.bokeh_chart(locations(df))
#UPS parameters input
    ups2 = pd.read_csv('ups2_dynamicStd.csv', sep=',')
    volume = st.number_input('How much volume - in µl - have you used to resuspend your UPS?', 0.0, 50.0, 21.2)
    amount = st.number_input('How much standard volume - in µl - have you spiked in your samples?', 0.0, volume, 6.0)
#UPS calculations based on the spiked amount of standards
    ups2, data, coef, inter = alpaca.ups(df, ups2, volume, amount)

    appended = list()
    #This new loop solves the problem of variable amount of experimental conditions.
    #It creates a dataframe (clean_data) based on the condition list.
    for x in condition:
        cond = alpaca.condition_format(df, x, rep)
        appended.append(cond)
    clean_data = pd.concat(appended)

    clean_data['Amount (fmol)'] = 2 ** ((clean_data[['Mean']] - inter) / coef)
    st.markdown('**Data preview with fmol amount of proteins based on UPS2 standard curve**')
    st.write(clean_data)

with expander2:
    st.header("Enrichment")
    enrichments = st.number_input('How many Enrichment conditions have your experiment?', 0, 10, 2)
    enrichment_t = list()
    en = 1
    for en in range(int(enrichments)+1):
        enrich_name = 'Enrichment_' + str(en+1)
        enrichment_t.append(enrich_name)
        en += 1

    if enrichments != 0 :
        instructions_place = st.empty()
        uploaded_std = st.file_uploader('Import your Enrichment Standards file',
                                        type = ['txt','csv','xlsx'],
                                    help='It takes a file for the calculation of the Enrichment factor.'
                                         'This table requires UNIPROT Accession, MW (kDa) and Mix concentration (µg/µl) '
                                         'as columns for the proper Enrichment calculation.')

        if uploaded_std is None:
            instructions_place.image(alpaca_enrichment)
            st.stop()

        std = alpaca.importer(uploaded_std)
        st.markdown('**Enrichment Standards**')
        stan = st.write(std)

        c1, c2, c3, c4 = st.columns(4)
        adder = 0
        enrichment_type_dict = dict()
        prep = list()
        for c in range(enrichments):
            enrich = c1.multiselect(enrichment_t[c], condition, condition[adder])
            dil = c2.number_input('Dilution', 1, 1000, 10 + adder)
            vol = c3.number_input('Spiked standard volume (µl)', 0.0, 1000.0, 8.5 + adder)
            sampl_v = c4.number_input('Sample volume (ml)', 0.0, 100.0, 45.0 + adder)
            prep.append(dil)
            prep.append(vol)
            prep.append(sampl_v)
            enrichment_type_dict[enrichment_t[c]] = enrich
            adder += 1

        data_std = alpaca.enricher(std, clean_data, prep, enrichments, n, condition)
        data_std = alpaca.enrichment_filter(data_std, enrichment_type_dict)
        data_std['Enrichment'] = data_std['Amount (fmol) in sample'] / data_std['Amount (fmol)']
        st.markdown('**Enrichment Standards in your data**')
        st.write('Data Dimension: ' + str(data_std.shape[0]) + ' rows and ' + str(data_std.shape[1]) + ' columns.')
        st.write(data_std)
        enrichment_factors = alpaca.enrichment_factors(data_std, enrichment_type_dict)
    elif enrichments == 0:
        columna1, columna2 = st.columns(2)
        adder = 0
        prep = list()
        for c in range(len(condition)):
            enrich = columna1.multiselect('Experimental conditions', condition, condition[adder])
            dil = 1
            vol = 1
            sampl_v = columna2.number_input('Sample volume (ml)', 0.0, 100.0, 45.0 + adder)
            prep.append(dil)
            prep.append(vol)
            prep.append(sampl_v)
            adder += 1
        enrichment_type_dict = dict()
        enrichment_type_dict['Non_Enrichment'] = condition
        enrichment_factors = dict()
        enrichment_factors['Non_Enrichment'] = 1

with expander3:

    st.header('Absolute Protein Quantification')
    st.balloons()

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

    dummy = alpaca.getEnriched(clean_data, enrichment_type_dict)
    dummy = alpaca.getMolecules(dummy, enrichment_factors)
    dummy = alpaca.getMolecules_cell(dummy, enrichment_type_dict, prep, count_dict)

    @st.cache
    def convert_df(df):
        # IMPORTANT: Cache the conversion to prevent computation on every rerun
        return df.to_csv().encode('utf-8')

    csv = convert_df(dummy)

    c1.download_button(
        label="Download data as CSV",
        data=csv,
        file_name='absolute_quantification.csv',
        mime='text/csv',
    )
    st.write('Data Dimension: ' + str(dummy.shape[0]) + ' rows and ' + str(dummy.shape[1]) + ' columns.')
    st.markdown('**Here comes the magic table**')
    AgGrid(dummy)
    st.stop()

"""
with expander4:
    st.header("Data visualisation")

    choice = st.selectbox('What would you like to check?',
                          ('ID proteins', 'Cellular locations', 'Cellular locations (Normalized)'))
    if choice == 'ID proteins':
        pepe = clean_data.groupby('Condition')[rep].count()
        pepe = pepe.reset_index()
        st.write(pepe)

    elif choice == 'Cellular locations':
        gp4 = pd.read_excel("201118_GP4_SCL_complete.xlsx")
        gp4[['sp', 'Accession', 'Entry Name']] = gp4['Protein_ID'].str.split('|', expand=True)
        gp4 = gp4.drop(columns='sp')

        loc = clean_data.merge(gp4, on='Accession', how='left')
        loc = loc[['Accession', 'T: Gene names', 'Predicted_SCL', 'Condition']]

        base = alt.Chart(loc).mark_bar().encode(
            x=('count(Predicted_SCL):Q'),
            y=('Predicted_SCL:N'),
            color=('Predicted_SCL:N'),
            row=('Condition:N')
        ).properties(
            width=600,
            height=200
        )
        st.altair_chart(base)

    elif choice == 'Cellular locations (Normalized)':
        gp4 = pd.read_excel("201118_GP4_SCL_complete.xlsx")
        gp4[['sp', 'Accession', 'Entry Name']] = gp4['Protein_ID'].str.split('|', expand=True)
        gp4 = gp4.drop(columns='sp')

        loc = clean_data.merge(gp4, on='Accession', how='left')
        loc = loc[['Accession', 'T: Gene names', 'Predicted_SCL', 'Condition']]

        base_n = alt.Chart(loc).mark_bar().encode(
            x=alt.X('count(Predicted_SCL):Q', stack='normalize'),
            y=('Predicted_SCL:N'),
            color=('Predicted_SCL:N'),
            row=('Condition:N')
        ).properties(
            width=600,
            height=200
        )
        st.altair_chart(base_n)




    #loc = gp4[gp4['Protein ID'].isin(df['Accession'])]

    #x = loc['Predicted_SCL'].astype(str).value_counts()
    #data = pd.Series(x).reset_index(name='value').rename(columns={'index':'location'})
    #data['angle'] = data['value']/data['value'].sum() * 2*pi
    #data['color'] = Category20c[len(x)]
    #z=110*(data['value']/data['value'].sum())
    #data['value']=z
    #p = figure(plot_height=350, title="", toolbar_location=None,
    #       tools="hover", tooltips="@location: @value{0.2f} %", x_range=(-.5, .5))
    #p.annular_wedge(x=0, y=1,  inner_radius=0.15, outer_radius=0.25, direction="anticlock",
    #            start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
    #    line_color="white", fill_color='color', legend='location', source=data)
    #p.axis.axis_label=None
    #p.axis.visible=False
    #p.grid.grid_line_color = None
    #st.bokeh_chart(p, use_container_width=False)
    #export_png(p, filename="location_plot.png")
    #with open("location_plot.png", "rb") as file:
    #    btn = st.download_button(
    #        label="Download graph",
    #        data=file,
    #        file_name="location_plot.png",
    #        mime="image/png"
    #    )
"""
