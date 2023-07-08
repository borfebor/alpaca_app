#needed packages
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image

from src.alpaca_new import alpaca
from src.viz import Viz

image = Image.open('ALPACA_LOGO2.png')
tab_logo = Image.open('tab_logo.png')

st.set_page_config(
     page_title="alpaca proteomics",
     page_icon=tab_logo,
     layout="wide",
     initial_sidebar_state="expanded"
)
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
    
    example_data = st.sidebar.selectbox('Example datasets', [None, 'Standard protocol', 'Enriched protocol'])

if example_data != None:   
    
    paper_dict = {'Standard protocol':'Cytosol_example', 
                  'Enriched protocol':'Enriched_example'}
    
    title = f'Working with an example dataset from {example_data}'
    top_bar.title(title)
    uploaded_file = f'Datasets/{paper_dict[example_data]}.txt'

if uploaded_file is None:
    st.image(instructions)
    st.stop()

if uploaded_file is not None:

    # Can be used wherever a "file-like" object is accepted:
    df = alpaca.eats(uploaded_file)
     
    st.sidebar.header('Data preprocessing')
    cleaning = st.sidebar.checkbox('Data cleaning', value=True)
    formatting = st.sidebar.checkbox('Data formatting', value=True)
    
    lfq_options = ['iBAQ', 'LFQ intensity', 'Top3', 'Intensity', 'MS/MS count']
    
    lfq_menu = list(dict.fromkeys([item for item in lfq_options for col in df.select_dtypes(include=np.number).columns if item in col]))
    
    lfq_method = st.sidebar.selectbox('Label-Free Quantification method', lfq_menu, 0)
    
    normalization = st.sidebar.selectbox('Intensity normalization', [False, 'Relative', 'Median', 'Quantile'], 0)
        
    df, condition, lfq_method = alpaca.spits(df, lfq_method=lfq_method, lfq_columns=lfq_menu,
                                                cleaning=cleaning, normalization=normalization,
                                                formatting=formatting, identifier=None)
    
    if formatting == False:
     
        export_not_formated_results = df.to_csv(sep='\t').encode('utf-8')
        
        st.dataframe(df, use_container_width=True)   
        
        st.download_button(
              label="Download clean data",
              data=export_not_formated_results,
              file_name='alpaca_clean_data.txt',
              mime='text/csv',
              help='Here you can download your data',
              use_container_width=True,
          )
        
        st.stop()
    
    st.sidebar.header('Quantification standards')
    
    standards = alpaca.eats('UPS2.txt')  
    
    used_std = st.sidebar.selectbox('Used quantification standards', 
                     ['UPS1', 'UPS2', 'UPS3', 'Custom'], 1)
    
    if used_std == 'Custom':
        
        custom_std = st.sidebar.file_uploader('Upload your custom standards')
        if custom_std == None:
          st.markdown("""
               ## Welcome to the Custom Standards assistant 🤖

               You can upload the used standards at the bottom of the sidebar, an uploader just poped up there. Just a few things to consider:

               1. Your standards should be included in the database used for the MaxQuant search, otherwise the standards will be missing in the dataset.
               2. Your standards list should contain the **same Accession** as the one they had in the database used for MaxQuant searches
               3. Files can be uploaded either as tab-delimitted (TXT), comma-delimitted (CSV) or excel (XLSX). I just need you to use the following format for preparing the standards list:

               |    Accession   |fmol                         	 |MW (Da)                         |
               |----------------|-------------------------------|-------------------|
               |P00689			 |5000     						 |20900       |
               |PX0709          |500           				 |49670            |
               |QZ06T9          |50 |18093|
               > **Amount column:** it should contain "fmol" on the title (e.g. Amount_fmol)
               
               > **Molecular weight:** It should be in Da and contain in the name either MW or Da for a proper identification (e.g. Molecular wight (Da)).
               """)
        else:
          standards = alpaca.eats(custom_std)  
        
    st.sidebar.subheader('Quantification standards parameters')
    
    samples = list(df.Sample.unique())
    replicate = None
        
    col1, col2 = st.sidebar.columns(2)  
    concentration = col1.number_input('Quantification Standards concentration (µg/µl)', 0.0, 21.2, 0.5)
    volume = 10.6 / concentration
    amount = col2.number_input('Added volume in sample (µL)', 0.0, volume, 6.0)
    
    samples_with_ups = st.sidebar.checkbox('All samples contain quantification standards', True) 
    
    standards_clean = standards.copy()
    
    if samples_with_ups == False:
        
        try:
            suggested = alpaca.pca(df, standards_clean, lfq_method)
            clustered, suggested = alpaca.KMeans(suggested)
        except:
            suggested = []
        
        replicate = st.sidebar.multiselect('Select where can I find your UPS spiked', 
                                   samples, suggested)
        
    selected = st.sidebar.multiselect('Is there any outlier to remove?',
                                list(standards['Accession'].unique()))
            
    standards_clean = standards[~standards['Accession'].isin(selected)] 
        
    experimenter = st.expander('Experimental set-up', expanded=False)
    
    enricher = st.expander('Enrichment standards', expanded=False)
        
    visualizer = st.expander('Data visualization', expanded=False)
    
    st.title('Your data')
    
    listed_cond = ', '.join(col for col in condition)
    st.write(f'Your data contains {len(condition)} experimental conditions ({listed_cond})')
    
    results = st.empty()
    
    data_frame = st.empty()
    
    try:                                                       
        data, ups2, coef, inter, R2 = alpaca.census(df, standards_clean, concentration=concentration, 
                                                           in_sample=amount, lfq_col=lfq_method, filter_col='Sample',
                                                           added_samples=replicate, valid_values=2,
                                                           save=False)   
    except:
        st.warning("Couldn't find the quantification standards")
        st.markdown("""
               ## Hi, I'm your Quantification Standards assistant 🤖

               Check that on the bottom part of the sidebar if the right standards are selected.
               Otherwise, you can upload your custom standards if you have not use any of the Universal Proteomics Standards mixes.
               
               Just a few things to consider:

               1. Your standards should be included in the database used for the MaxQuant search, otherwise the standards will be missing in the dataset.
               2. Your standards list should contain the **same Accession** as the one they had in the database used for MaxQuant searches
               3. Files can be uploaded either as tab-delimitted (TXT), comma-delimitted (CSV) or excel (XLSX). I just need you to use the following format for preparing the standards list:

               |    Accession   |fmol                         	 |MW (Da)                         |
               |----------------|-------------------------------|-------------------|
               |P00689			 |5000     						 |20900       |
               |PX0709          |500           				 |49670            |
               |QZ06T9          |50 |18093|
               > **Amount column:** it should contain "fmol" on the title (e.g. Amount_fmol)
               
               > **Molecular weight:** It should be in Da and contain in the name either MW or Da for a proper identification (e.g. Molecular wight (Da)).
               """)
        st.stop()
    
    with experimenter:
        
          st.markdown("""
            #### Welcome to the experimental details panel! 👋
            Here you will be able to translate the labwork into data analysis.
            """,
            )
          
          tab1, tab2, tab3 = st.tabs(["📥 Upload your experimental details", 
                                      "🔬 Check your experimental parameters", 
                                      "📍 Check your enrichment standards"])
          
          with tab1:
              
              st.markdown("""
                          Here you can upload your sample parameters.
                          
                          1. **Experimental parameters:** It should contain the details of your sample preparation. 
                          File name should contain **'param'**.
                          2. **Enrichment standards:** Just needed if there has been an enrichment step. 
                          It should contain the details of the used protein standards. File name should contain **'standard'**.
                          
                          Examples can be downloaded or previewed by clicking the buttons below
                         """)
                         
              table = st.empty()
              
              example = st.empty()
              
              params_temp, params_show, std_temp, std_show = st.columns([1, 1, 1, 1])

              params_example = alpaca.generate_example_params(df).reset_index(drop=True)
              
              if params_show.button('Show params example'):

                  example.dataframe(params_example.set_index('Condition'), use_container_width=True)
                  
              enrichment_example = alpaca.create_random_df(data).reset_index(drop=True)
              
              if std_show.button('Show standards example'):
                  
                  example.dataframe(enrichment_example, use_container_width=True)
              
              export_params = params_example.to_csv(sep='\t').encode('utf-8')
          
              params_temp.download_button(
                    label="Get params template",
                    data=export_params,
                    file_name='alpaca_params.txt',
                    mime='text/csv',
                    help='I will help you',
                    use_container_width=True,
                )
              
              export_stds = enrichment_example.to_csv(sep='\t').encode('utf-8')
          
              std_temp.download_button(
                    label="Get standards template",
                    data=export_stds,
                    file_name='alpaca_standards.txt',
                    mime='text/csv',
                    help='I will help you',
                    use_container_width=True,
                )
                           
              experiment = table.file_uploader('Import your params and/or standards file',
                                                  type=['csv','txt','xlsx'],
                                                   help="""Click here to upload your experimental parameters 
                                                   and/or enrichment standards. Accepted file formats are Excel (.xlsx), 
                                                   CSV (.csv), and TXT (.txt). If uploading enrichment standards, 
                                                   please name the file "standards".
                                                      """,
                                            accept_multiple_files=True)      
          
          
              for item in experiment:
                  
                  if "PARAM" in item.name.upper():
                      
                      sample_conditions = alpaca.eats_better(item)
                      sample_conditions = alpaca.matchmaker(sample_conditions, params_example)
                      sample_conditions = sample_conditions.set_index(sample_conditions.columns[0])
                      
                      #agree = agree_placeholder.checkbox('Use parameters for calculations', True)
                      
              for item in experiment:
                      
                  if "STANDARD" in item.name.upper():
                      
                      enrichment_std = alpaca.eats(item)
                      enrichment_std = alpaca.matchmaker(enrichment_std, enrichment_example)
                      
          with tab2:
                
              placer_tab3 = st.empty()
              
              if 'sample_conditions' in globals():
                    
                    placer_tab3.dataframe(sample_conditions)
                
              else:
                    placer_tab3.markdown("""
                                Someone is trying to run too fast.
                                Please, upload your experimental parameters file first.
                               """)
                               
          with tab3:
              
              placer_tab4 = st.empty()
              if 'enrichment_std' in globals():
                  if 'sample_conditions' in globals():
                      
                      e_test, sample_conditions = alpaca.gathers(data, enrichment_std, sample_conditions, lfq_method=lfq_method) 
                      placer_tab4.dataframe(e_test, use_container_width=True)
                  
              else:
                      placer_tab4.markdown("""In case you have performed an enrichment step in your sample preparation,
                                 do not forget to upload the enrichment standards. 
                                 Please, find an example below on how the input should look like.
                                 """)
                                       
       
    if 'sample_conditions' in globals():
        
        data = alpaca.wooler(df, sample_conditions)


    param_1 = 'Sample'
    param_2 = 'iBAQ'
    param_3 = 'Condition'
     
    with visualizer:
        
        settings, plot = st.columns([1,2])
        
        plots = ['Intensities', 'Quantified proteins', 'Calibration curve',
         'PCA','Distribution plot', 'Heatmap']
            
        if 'e_test' in globals():
            plots.append('Enrichment')
            
        viz_type = settings.selectbox('What would you like to inspect?', 
                                plots, 0)
        
        adjusment = settings.checkbox('Adjust plot to the container', True)
        
        categories = [col for col in df.select_dtypes('object').columns]   
        numbers = [col for col in df.select_dtypes('number').columns]
        
        sindex = categories.index('Sample')
        con_index = categories.index('Condition')
        lfq_index = numbers.index(lfq_method)
        
        if viz_type == 'Intensities':
                               
                categorical = settings.selectbox('X axis category', 
                                                 categories, sindex)
                numerical = settings.selectbox('Y axis values', 
                                               numbers, lfq_index)
                
                color = settings.selectbox('Pick a category', 
                                                 categories, con_index)
             
                box = Viz.boxplot_2(df, categorical, numerical, color)
                #chart  = Viz.boxplot(df, categorical, numerical, color)
                plot.plotly_chart(box)
                
        elif viz_type == 'Quantified proteins':
            
            grouper = settings.multiselect('Group your data', 
                                             categories, ['Sample', 'Condition'])
            
            categorical = settings.selectbox('X axis category', 
                                             grouper, 0)
            
            color = settings.selectbox('Pick a category', 
                                             grouper, 1)       
            
            chart  = Viz.identifications(df, categorical, color, grouper)
                
        elif viz_type == 'Calibration curve':
                       
            settings.write(ups2)
            
            chart  = Viz.Regression(ups2, 'log2_Amount_fmol', lfq_method, R2)
            
        elif viz_type == 'PCA':
            
                index = categories.index('Condition')
            
                sorter = settings.selectbox('Filtering category', categories, index)
                sorting = settings.multiselect('Exclude items', df[sorter].unique()) 
                
                to_pca = df[~df[sorter].isin(sorting)]
                
                pcanalysis, components, variance = alpaca.pcomponent(to_pca, lfq_method, 
                                                                     conditions=condition)
                
                x = settings.selectbox('X axis component', 
                                                 components, 0)
                
                y = settings.selectbox('Y axis component', 
                                                 components, 1)
                
                chart  = Viz.scatter(pcanalysis, x=x, y=y, variance=variance, columns=components)
                
        elif viz_type == 'Distribution plot':
            
                chart  = Viz.displot(df, lfq_method)
                
        elif viz_type == 'Enrichment':
            
                grouper = settings.multiselect('Group your data', 
                                                 categories, ['Condition', 'Replicate'])
                
                categorical = settings.selectbox('X axis category', 
                                                 grouper, 0)
                
                color = settings.selectbox('Pick a category', 
                                                 grouper, 0)  
                
                e_group = e_test.groupby(grouper)['Enrichment'].median().reset_index()
            
                chart = Viz.boxplot(e_group, categorical, 'Enrichment', color)
                
        elif viz_type == 'Heatmap':
                     
            z_score = settings.checkbox('Z score normalization', False)
            
            numerical = settings.selectbox('Value to plot', 
                                           numbers, lfq_index)
            
            color_scheme = settings.selectbox('Colour scheme', 
                                           ['redblue', 'blueorange', 'spectral', 'tealblues', 'purples', 
                                            'lightorange', 'lighttealblue'], 0)
            
            chart = Viz.heatmap(data, 'Sample','Protein',  numerical, z_score, color_scheme)

        plot.altair_chart(chart, use_container_width=adjusment) 
    
    with results:
                    
        c1, c2, c3 = results.columns([1, 2, 3])
        
        if c1.checkbox('Pivot table'):
        
            data = alpaca.pivoter(data)
            
        select_cols = [col for col in data.select_dtypes(exclude=np.number)]
            
        look_on = c2.selectbox('Search based on', select_cols)
     
        query = c3.multiselect('Search term (e.g., TrxA, Erv1, Cox1)', data[look_on].unique(),
                               default=[])
        
        if query != []:
            
            data = data[data[look_on].isin(query)]
        
        export_results = data.to_csv(sep='\t').encode('utf-8')
        
        data_frame.dataframe(data, use_container_width=True)   
        
        c1.download_button(
              label="Download results",
              data=export_results,
              file_name='alpaca_results.txt',
              mime='text/csv',
              help='Here you can download your results',
              use_container_width=True,
          )

