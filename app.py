#needed packages
import numpy as np
import streamlit as st
from PIL import Image
import pandas as pd

from src.alpaca_st import alpaca
from src.viz import Viz
from src.spits import detective, Normalization, Imputation, tools
from src.census import advisor

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
                                        type=['csv','txt','xlsx', 'tsv'],
                                         help='You can upload directly your ProteinGroup.txt from MaxQuant output '
                                              'and clean it using the button below, or pre-clean the data with Perseus. '
                                              'Different intensity-based quantification methods can be selected in the app.')

example_data = None

if uploaded_file == None:
    
    example_data = st.sidebar.selectbox('Example datasets', [None, 'Ferrero-Bordera et al. 2024', 'Antelo-Varela et al. 2019'])
    st.sidebar.markdown("""
                ### [Tutorial and documentation](https://github.com/borfebor/alpaca_app)
                
                Developed by [Microbial Proteomics (Uni Greifswald)](https://microbialproteomics.uni-greifswald.de/en/)
                
                """)
                
if example_data != None:   
    
    paper_dict = { 
                  'Ferrero-Bordera et al. 2024':{'id':'Enriched_example','source':'Ferrero-Bordera et al. 2024. Microbiology Spectrum','link':'https://doi.org/10.1128/spectrum.02616-23'},
                  'Antelo-Varela et al. 2019':{'id':'Membrane_example','source':'Antelo-Varela et al. 2019. Anal. Chem.','link':'https://doi.org/10.1021/acs.analchem.9b02869'}
                  }
    
    top_bar.markdown(f"Working with an example dataset from [{paper_dict[example_data]['source']}]({paper_dict[example_data]['link']})")
    uploaded_file = f"Datasets/{paper_dict[example_data]['id']}.txt"

if uploaded_file is None:
    st.image(instructions)
    st.stop()

if uploaded_file is not None:

    # Can be used wherever a "file-like" object is accepted:
    df = alpaca.eats(uploaded_file)

    aid =  st.empty()
    
    st.sidebar.header('Quantification standards')
    
    standards = alpaca.eats('UPS2.txt')
    
    used_std = st.sidebar.selectbox('Used quantification standards', 
                     ['UPS2', 'Custom'], 0)
    
    if used_std == 'Custom':
        
        custom_std = st.sidebar.file_uploader('Upload your custom standards')
        if custom_std == None:
          aid.markdown("""
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
               
               > **Molecular weight:** It should be in Da and contain in the name either MW or Da for a proper identification (e.g. Molecular weight (Da)).
               """)
               
        else:
          standards = alpaca.eats(custom_std)  
        

     
    st.sidebar.header('Data preprocessing')
    cleaning = st.sidebar.checkbox('Data cleaning', value=True)
    
    id_col, is_pivot, it = detective.alpacaHolmes(df)
    it_raw = it.copy()
    
    formatting = st.sidebar.checkbox('Data formatting', value=is_pivot)
    
    lfq_options = list(it.keys())
    
    id_cols = list(dict.fromkeys([col for col in df.select_dtypes(exclude=[np.number, 'bool']).columns]))
    lfq_menu = list(dict.fromkeys([item for item in lfq_options for col in df.select_dtypes(include=np.number).columns if item in col]))
    
    id_col = st.sidebar.selectbox('Accession column', id_cols, id_cols.index(id_col))
    lfq_method = st.sidebar.selectbox('Label-Free Quantification method', lfq_menu, 0)
    
    replicate_dict, conditions, it = detective.alpacaWatson(df, it, id_col, lfq_method)
    
    norm_dict = {'None':None, 'Relative': Normalization.Relative, 'Median': Normalization.Median, 'Quantile': Normalization.Quantile}
    
    normalization = st.sidebar.selectbox('Intensity normalization', list(norm_dict.keys()), 0)
    
    the_list = [len(v) for k, v in  tools.invert_dict(replicate_dict).items()]
    max_values = np.max(the_list)
    
    value_filter = st.sidebar.slider('Valid values per condition', min_value=0, max_value=max_values, value=2) / max_values
    
    imputation_methods = {
        'None': ('', {''}),
        'Lowest of Detection (LOD)': (Imputation.impute_lod, {'lod': 0.01}),
        'Normal Distribution (ND)': (Imputation.impute_nd, {'lod': 0.01}),# 'mean': mean, 'std': std}),
        'k-Nearest Neighbors (kNN)': (Imputation.impute_knn, {'n_neighbors': 5}),
        'Local Least Squares (LLS)': (Imputation.impute_lls, {'max_iter': 10}),
        #'Random Forest (RF)': (Imputation.impute_rf, {'max_iter': 10}),
        'Singular Value Decomposition (SVD)': (Imputation.impute_svd, {'n_components': 2}),
        #'Bayesian Principal Component Analysis (BPCA)': (Imputation.impute_bpca, {'n_components': 2, 'max_iter': 100})
        }
    
    impute = st.sidebar.selectbox('Imputation method', list(imputation_methods.keys()), 0)

    imp_method = imputation_methods[impute][0]
    
    st.sidebar.subheader('Quantification standards parameters')
        
    col1, col2 = st.sidebar.columns(2)  
    concentration = col1.number_input('Quantification Standards concentration (µg/µl)', 0.0, 21.2, 0.5)
    volume = 10.6 / concentration
    amount = col2.number_input('Added volume in sample (µL)', 0.0, volume, 6.0)
    
    samples_with_ups = st.sidebar.checkbox('All samples contain quantification standards', True) 
    
    standards_clean = standards.copy()
    
    raw_df = df.copy()
    df = alpaca.spits(df, id_col, lfq_method, replicate_dict,
                      cleaning=cleaning, formatting=formatting, normalization=norm_dict[normalization],
                      valid_values=value_filter, imputation=imp_method,)
    
    samples = list(replicate_dict.values())
   
    replicate = None
    added = None
    
    if samples_with_ups == False:
        
        try:
            suggested = advisor.pca(df, standards_clean, lfq_method)
            clustered, suggested = advisor.KMeans(suggested)
        except:
            suggested = []
        
        suggested = [replicate_dict[i] for i in suggested]
        added = st.sidebar.multiselect('Select where can I find your UPS spiked', 
                                   samples, suggested)
        
        inverted = {v: k for k, v in replicate_dict.items()}
        
        replicate = [inverted[i] for i in added]

        
    selected = st.sidebar.multiselect('Is there any outlier to remove?',
                                list(standards['Accession'].unique()))
            
    standards_clean = standards[~standards['Accession'].isin(selected)] 
    
    if formatting == False:
     
        export_not_formated_results = df.to_csv(sep='\t').encode('utf-8')
        
        st.title('Your data')
        
        listed_cond = ', '.join(col for col in conditions)
        st.write(f'Your data contains {df.Accession.nunique()} proteins from {len(conditions)} experimental conditions ({listed_cond})')
        
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
        
    experimenter = st.expander('Experimental set-up', expanded=False)
    
    enricher = st.expander('Enrichment standards', expanded=False)
        
    visualizer = st.expander('Data visualization', expanded=True)
    
    st.title('Your data')
    
    listed_cond = ', '.join(col for col in conditions)
    st.write(f'Your data contains {df.Accession.nunique()} proteins from {len(conditions)} experimental conditions ({listed_cond})')
    
    message = st.empty()
    
    results = st.empty()
    
    data_frame = st.empty()
    
    try:                                                       
        data, ups2, coef, inter, R2 = alpaca.census(df, standards_clean, concentration=concentration, 
                                                           in_sample=amount, lfq_col=lfq_method, filter_col='Sample',
                                                           added_samples=replicate, valid_values=2,
                                                           save=False)  
        if R2 < 0.8:
            message.warning(f"Low R² value ({R2:.2f}), indicating a poor fit for the regression.")
    except:
        aid.warning("Couldn't find the quantification standards")
        aid.markdown("""
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
        with results:
        
                      
                      export_results = df.to_csv(sep='\t').encode('utf-8')
                      
                      data_frame.dataframe(df, use_container_width=True)   
                      
                      st.download_button(
                            label="Download results",
                            data=export_results,
                            file_name='alpaca_results.txt',
                            mime='text/csv',
                            help='Here you can download your results',
                            use_container_width=True,
                        )
               
    try:
        exist = ups2
    except:
        exist = None
    
    if exist is not None:
           
        std_test = ups2[['Accession', 'log2_Amount_fmol']].copy()
        
        norm_advise = ''
        res_hm, lfq_advise, norm_advise = advisor.mini_spits(raw_df, std_test, it, id_col, norm_dict=norm_dict,
                           added_samples=added)
        
        if norm_advise != '':
            aid.info(f'Based on your data, {norm_advise}-normalized {lfq_advise} is recommended for the quantification.')
        else:
            aid.info(f'Based on your data, {lfq_advise} without normalization is recommended for the quantification.')

    
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
                  
              enrichment_example = alpaca.create_random_df(df).reset_index(drop=True)
              
              if std_show.button('Show standards example'):
                  
                  example.dataframe(enrichment_example, use_container_width=True)
              
              export_params = params_example.to_csv(sep='\t', index=False).encode('utf-8')
          
              params_temp.download_button(
                    label="Get params template",
                    data=export_params,
                    file_name='alpaca_params.txt',
                    mime='text/csv',
                    help='I will help you',
                    use_container_width=True,
                )
              
              export_stds = enrichment_example.to_csv(sep='\t', index=False).encode('utf-8')
          
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
                      
                      sample_conditions = alpaca.eats(item)
                      sample_conditions = alpaca.matchmaker(sample_conditions, params_example, 85)
                      sample_conditions = sample_conditions#.set_index(sample_conditions.columns[0])
                      
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
        
        data = alpaca.wool(df, sample_conditions)


    param_1 = 'Sample'
    param_2 = 'iBAQ'
    param_3 = 'Condition'
     
    with visualizer:
        
        settings, plot = st.columns([1,2])
        
        plots = []
                
        if 'res_hm' in globals():
            plots.append('Fitting advisor')
            
        plots = plots + ['Intensities', 'Quantified proteins', 'Calibration curve',
         'PCA','Distribution plot', 'Heatmap']
        
        if 'e_test' in globals():
            plots.append('Enrichment')
            
        viz_type = settings.selectbox('What would you like to inspect?', 
                                plots, plots.index('Quantified proteins'))
        
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
             
                chart = Viz.boxplot(df, categorical, numerical, color)
                
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
                                                                     conditions=conditions)
                
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

            hm_palettes = ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
                                            'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
                                            'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
                                            'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
                                            'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
                                            'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
                                            'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                                            'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
                                            'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
                                            'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
                                            'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
                                            'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
                                             'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
                                            'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr','ylorrd']
            
            color_scheme = settings.selectbox('Colour scheme', 
                                         hm_palettes , hm_palettes.index('rdbu'))
            
            chart = Viz.heatmap(data, 'Sample','Protein',  numerical, z_score, color_scheme)
            
        elif viz_type == 'Fitting advisor':
                     
            z_score = False
            res_hm['score'] = res_hm['score'].apply(lambda x: np.round(x, 3))
            hm_palettes = ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
                                            'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
                                            'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
                                            'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
                                            'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
                                            'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
                                            'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
                                            'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
                                            'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
                                            'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
                                            'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
                                            'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
                                             'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
                                            'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr','ylorrd']
            
            color_scheme = settings.selectbox('Colour scheme', 
                                         hm_palettes , hm_palettes.index('blues'))
            
            chart = Viz.heatmap(res_hm, 'method', 'normalization',  'score', z_score, color_scheme)

        try: 
             plot.plotly_chart(chart, theme=None, use_container_width=adjusment) 
        except:
             plot.error(chart, icon="🚨")
    
    if 'data' not in globals():
         st.stop()
     
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

     
