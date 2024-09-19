import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import re
from src.spits import Normalization, tools, clean, Imputation
from src.census import Quantification
from src.gathers import gathering
from src.wool import yarnScissors

class alpaca:
    
    def eats(file):
        """
        Loads a file (either in CSV, TSV, or XLSX format) into a Pandas DataFrame.
        
        Parameters:
        file : str or file-like object
            The file can be provided as a string representing the file path (in which case it's assumed to be on disk)
            or as a file-like object (e.g., uploaded file in Streamlit).
        
        Returns:
        df : pandas.DataFrame or None
            The file content as a DataFrame if it is of a supported type (CSV, TSV, XLSX), otherwise returns None.
        """
        
        if isinstance(file, str):
            file_name = file
        else:
            file_name = file.name

        if 'TXT' in file_name.upper() or 'TSV' in file_name.upper():
            df = pd.read_csv(file, sep='\t')  # Tab-separated values
        elif 'CSV' in file_name.upper():
            df = pd.read_csv(file, sep=',')  # Comma-separated values
        elif 'XLSX' in file_name.upper():
            df = pd.read_excel(file)  # Excel file
        else:
            st.warning('Not compatible format')
            return None  # Return None explicitly if format is unsupported
    
        return df
       
    def spits(df, id_col, lfq_method, replicate_dict, cleaning=True, formatting='auto', 
             transformation=np.log2, normalization=None, valid_values=0.7, imputation='', **imp_kwargs):
        """
        Processes a DataFrame for quantitative analysis, performing tasks such as cleaning, 
        transformation, imputation, normalization, and reformatting based on provided parameters.
    
        Parameters:
        df : pandas.DataFrame
            Input dataframe containing raw data for processing.
        id_col : str
            Column name in `df` representing unique identifiers for the dataset (e.g., 'Accession').
        lfq_method : str
            Column name or method representing LFQ (Label-Free Quantification) values in the data.
        replicate_dict : dict
            A dictionary mapping sample names to condition and replicate information.
        cleaning : bool, optional, default=True
            Whether to clean the data by removing contaminants and decoys.
        formatting : str or bool, optional, default='auto'
            Controls the format of the output dataframe ('auto', True, or False).
        transformation : callable, optional, default=np.log2
            A transformation function to apply to the data (e.g., log2).
        normalization : callable or None, optional, default=None
            A normalization function to apply to the data. If None, no normalization is applied.
        valid_values : float, optional, default=0.7
            Proportion of valid values required for each row to be retained.
        imputation : str, optional, default=''
            Method for data imputation. If empty, no imputation is performed.
        **imp_kwargs : dict
            Additional keyword arguments passed to the imputation function.
        
        Returns:
        pandas.DataFrame
            Processed DataFrame ready for further analysis or visualization.
        """
        
        df.columns = df.columns.str.replace('.: ', '')
        samples = list(replicate_dict.keys())
        #conditions = list(tools.invert_dict(replicate_dict).keys())
    
        if 'Accession' not in df.columns:
            df.columns = [re.sub(id_col, 'Accession', id_col) 
                          if id_col == col else col for col in df.columns]
    
        # Removal of contaminants, and decoys
        potential_cols = ['identified by site', 'contaminant', 'Reverse']
        cont_key = [col for col in df.select_dtypes(exclude=np.number).columns for item in potential_cols if item in col]
    
        # Select columns with additional information 
        all_ids, priority = ['Accession'],  ['name', 'kDa']
        [all_ids.append(col) for col in df.columns for i in priority if i.upper() in col.upper()]
        [all_ids.append(col) for col in df.select_dtypes(exclude=np.number).columns if col not in all_ids]
        wanted_ids = st.sidebar.multiselect('Data identifiers of interest', all_ids , all_ids[:3])
        ids = [col for col in df.select_dtypes(exclude=np.number).columns if col in wanted_ids]
        
        if cleaning is True:
            cont_key = st.sidebar.multiselect('Items to remove', cont_key, cont_key)
            df = df[df[cont_key].isna().all(axis=1)].reset_index(drop=True)
            print(f"Items marked on {', '.join(cont_key)} have been removed from the dataset.")
    
        df = df[ids+samples].replace(0, np.nan)
        log_samples = [i for i in samples if 'LOG' not in i.upper()]
        df[log_samples] = df[log_samples].apply(lambda x: transformation(x)) 
        df[samples] = clean.filter_rows_by_missing_values(df[samples], replicate_dict, valid_values)
        #df = df.dropna(subset=samples, thresh=1)
        if imputation != '':
            df[samples] = Imputation.by_condition(df[samples], replicate_dict, imputation, **imp_kwargs)
        df = Normalization.normalizer(df, samples, normalization)
        
        if formatting == 'auto':
            formatting = tools.check_format(df, id_col)
    
        if formatting is True:
            df = df.melt(id_vars=ids, value_vars=samples, 
                        var_name='Sample', value_name=lfq_method).replace(0, np.nan)
            df[['Condition', 'Replicate']] = df['Sample'].replace(replicate_dict).str.split(';', expand=True)
            df = df.dropna(subset=lfq_method)

            if 'Gene names' in ids:
                df['Gene names'] = df['Gene names'].str[0].str.upper() + df['Gene names'].str[1:]
                df = df.rename(columns={'Gene names':'Protein'})
                
            print('Dataset formated for further analysis and visualisation.')   
            return df
        else:
            st.warning('Data is formated for human vision.\nThat could lead to errors or incompatibilities in further analysis using Alpaca pipeline.\nConsider formating your dataset if you see any anomally.')
            return df

        
    def census(df, standards, concentration=0.5, in_sample=6.0, lfq_col='iBAQ', ratio=1, 
               total_protein= 1, filter_col='Replicate', added_samples=None, valid_values=2, save=True):
        '''
        Performs protein quantification using regression between quantified protein intensities and 
        dynamic standards (e.g., UPS2).
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing clean data for quantified proteins.
        standards : str or file-like object (.csv, .txt, .xlsx)
            Path to a file containing UPS2 dynamic standards information or the standards file itself.
        concentration : float, optional, default=0.5
            Stock concentration of the standards in mg/mL.
        in_sample : float, optional, default=6.0
            Volume (in microliters) of standards added to each sample.
        lfq_col : str, optional, default='iBAQ'
            Column in `df` representing the label-free quantification (LFQ) values. Can be 'iBAQ', 'LFQ', or 'Intensity'.
        ratio : float, optional, default=1
            A multiplier for adjusting the calculated concentration of each protein.
        total_protein : float, optional, default=1
            Total protein concentration in the sample.
        filter_col : str, optional, default='Replicate'
            Column to filter data for specific replicates or conditions.
        added_samples : list or None, optional, default=None
            List of samples or conditions where standards were added. If None, assume standards were added to all samples.
        valid_values : int, optional, default=2
            Minimum number of valid (non-missing) values required for regression.
        save : bool, optional, default=True
            Whether to save regression plots.
    
        Returns
        -------
        df : pandas.DataFrame
            DataFrame with quantified proteins, adjusted for concentration.
        ups_red : pandas.DataFrame
            DataFrame containing measured standards in the sample.
        coef : float
            Regression slope (used for quantification).
        inter : float
            Regression intercept (used for quantification).
        R2 : float
            R-squared value representing the goodness-of-fit for the regression.
            
        '''
        try:
            ups2 = Quantification.abacus(standards, concentration, in_sample, total_protein=total_protein)
        except Exception as e:
            raise ValueError(f"Error processing standards: {e}")
        if lfq_col not in df.columns:
            raise ValueError(f"The specified lfq_col '{lfq_col}' does not exist in the dataframe.")
        if added_samples is None:
            added_samples = df[filter_col].unique().tolist()  # Apply to all replicates by default
        ups_red, coef, inter, R2 = Quantification.regression(df, ups2, lfq_col=lfq_col, filter_col=filter_col,
                                                     added_samples=added_samples, valid_values=valid_values) # Regression between intensities and standards
        #if R2 < 0.8:
        #    st.warning(f"Low R² value ({R2:.2f}), indicating a poor fit for the regression.")
        Quantification.abacusreg(ups_red, lfq_col=lfq_col, R2=R2, save=save) # Plots the Regression
        df = Quantification.moles(df, coef, inter, lfq_col=lfq_col, ratio = ratio)
    
        return df, ups_red, coef, inter, R2
    
    def gathers(df, enrichment_standards, preparation, plot=False, save_plot=False, lfq_method='iBAQ'):
        """
        Calculates enrichment factors for protein standards spiked into samples and optionally plots the results.
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the quantified protein data (e.g., with iBAQ or LFQ values).
        enrichment_standards : pandas.DataFrame
            DataFrame or file containing the standard information for calculating enrichment.
        preparation : pandas.DataFrame
            DataFrame containing information about the experimental preparation (e.g., conditions).
        plot : bool, optional
            Whether to plot the enrichment factors. Default is False.
        save_plot : bool, optional
            Whether to save the plot as a file. Default is False.
        lfq_method : str, optional
            Label-free quantification method used in the analysis (e.g., 'iBAQ', 'LFQ'). Default is 'iBAQ'.
    
        Returns
        -------
        e_test : pandas.DataFrame
            DataFrame containing calculated enrichment values for each sample.
        preparation : pandas.DataFrame
            Updated preparation DataFrame with enrichment factors added.
        """
            
        e_test = gathering.enrichment_calculator(df, enrichment_standards, preparation, lfq_method).dropna(subset='Enrichment')
        required_cols = ['Condition', 'Replicate', 'Enrichment']
        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"Missing required columns in DataFrame: {required_cols}")
        grouping = ['Condition', 'Replicate']
        col_grouper = [columns for columns in df.columns if columns in grouping]

        enrichments = e_test.groupby(col_grouper)['Enrichment'].median().reset_index()
        col_grouper.remove('Replicate')
        enrichment_factors = enrichments.groupby(col_grouper).apply(lambda x: pd.Series({
                                                                    'EnrichmentFactor':x["Enrichment"].median()})
            ).reset_index()

        preparation = preparation.merge(enrichment_factors, on='Condition', how='left')
            
        for index, row in enrichment_factors.iterrows():
            print(f'Enrichment factor for condition {row["Condition"]}: {row["EnrichmentFactor"]}')
     
        if plot == True:
            sns.catplot(data=enrichments, x='Condition', y='Enrichment', kind='box', width=0.5)
        if save_plot == True:
            plt.savefig('spiked_in_standards.svg', bbox_inches='tight', pad_inches=0.5)
        
        return e_test, preparation
   
    def preparator(std, clean, sample, lfq_method='iBAQ'):
    
        enriched_conditions = [condition for condition, details in sample['Added volume'].items() if details != 0]
    
        std_conditions = pd.DataFrame()
        for condition in enriched_conditions:
            std['Condition'] = pd.Series([condition]*std.shape[0])
            std_conditions = pd.concat([std_conditions, std])
    
        e_test = std_conditions.copy()
    
        e_test['Added_ng'] = np.nan
        e_test['Added_fmol'] = np.nan
    
        MW = [col for col in e_test.columns if 'Da' in col][0]
    
        for condition in enriched_conditions:
    
            e_test['Added_ng'] = np.where(e_test.Condition == condition, 
                                          e_test['Mix concentration (µg/µl)'] / sample['Dilution'][condition] * sample['Added volume'][condition], 
                                          e_test['Added_ng'])
    
            e_test['Added_fmol'] = np.where(e_test.Condition == condition, 
                                          e_test['Added_ng'] / e_test[MW] * 1000, 
                                            e_test['Added_fmol'])
    
            e_test['Added_fmol'] = np.where(e_test['Added_fmol'] == 0, 
                                          np.nan, e_test['Added_fmol'])
    
        e_reduction = ['Accession', 'Mol. weight [kDa]', 'Sample', lfq_method, 
                               'Condition', 'Replicate', 'fmol','Added_fmol', 'Subproteome']
    
        e_test = e_test[[col for col in e_test.columns if col in e_reduction]]
    
        e_found = e_test.merge(clean, on=['Accession', 'Condition'], how='left')
        e_found['Enrichment'] = e_found['Added_fmol'] / e_found['fmol']
        
        return e_found    
     
    def wool(df, preparation):
        """
        Applies enrichment factor corrections to a dataframe, converts amounts to molecules using Avogadro's number,
        and computes sample-specific and cell-specific molecule concentrations.
    
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing the quantified protein data (e.g., fmol values).
        preparation : pandas.DataFrame
            DataFrame containing experimental setup and enrichment information (e.g., Enrichment factors, sample volumes).
    
        Returns
        -------
        df : pandas.DataFrame
            Updated DataFrame with corrected fmol values and molecule counts.
        """
        enrichment_params = ['Enrichment', 'EnrichmentMode', 'ProteinSRM', 'fmolSRM', 'EnrichmentFactor']
        sample_params = ['SampleVolume', 'ProteinConcentration', 'AmountMS']
        cells_params = ['CellsPerML', 'TotalCultureVolume']
        
        required_columns = ['Condition', 'fmol']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"Required columns {required_columns} are missing in the input DataFrame")
        
        if 'EnrichmentFactor' in preparation.columns:
        
            for condition, values in preparation.set_index('Condition')[enrichment_params].fillna(1).iterrows():
                
                if values['EnrichmentMode'] == 'Amplification':
                    """
                    This calculation is made for samples which correspond to a higher fraction 
                    compared to the original proteome. E.g., Membrane
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                        df.fmol / values['EnrichmentFactor'], 
                                        df.fmol)
                    
                elif values['EnrichmentMode'] == 'Sampling': 
                    """
                    This calculation is made for samples which correspond to a smaller fraction
                    to the original proteome. E.g., Secretome
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                          df.fmol * values['EnrichmentFactor'], 
                                          df.fmol)
                    
            preparation = yarnScissors.correctionSRM(df, preparation)

            if "CorrectionSRM" in preparation.columns:
                
                for condition, values in preparation.set_index('Condition').fillna(1).iterrows():
                        
                            df['fmol'] = np.where(df.Condition == condition, 
                                                df.fmol * values['CorrectionSRM'], 
                                                df.fmol)
        
        df['Molecules'] = df['fmol'] * 6.023e8  # Avogadro's number fixed for fmol (-15)
        
        if all(item in preparation.columns.to_list() for item in sample_params):
            
            df['fmolSample'] = np.nan
            for condition, values in preparation.set_index('Condition')[sample_params].fillna(1).iterrows():
                
                total_protein = values['SampleVolume'] * values['ProteinConcentration'] #calculate the µg of protein in the sample
                MS_to_sample = total_protein / values['AmountMS']
                
                df['fmolSample'] = np.where(df['Condition'] == condition,
                                                    df['fmol'] * MS_to_sample, 
                                                    df['fmolSample'])
                
                df['Molecules'] = df['fmolSample'] * 6.023e8 # Avogadro's number fixed for fmol (-15)
            
        if all(item in preparation.columns.to_list() for item in cells_params):
               
            df['MoleculesPerCell'] = np.nan
            for condition, values in preparation.set_index('Condition')[cells_params].fillna(1).iterrows():
                
                cells = values['CellsPerML'] * values['TotalCultureVolume']
                
                df['MoleculesPerCell'] = np.where(df['Condition'] == condition,
                                                    df['Molecules'] / cells, 
                                                    df['MoleculesPerCell'])
        
        return df
    
    def scienttist(prep, enrichment_type_dict, subproteome_dict=None):
        preparation = dict()
        for prepa in enumerate(enrichment_type_dict):
            enriched = prep[prepa[0]][1] != 0
            if subproteome_dict != None:
                preparation[prepa[1]]= {'Enrichment': enriched,
                                 'Dilution': prep[prepa[0]][0],
                                 'Added_vol': prep[prepa[0]][1],
                                 'Sample_vol': prep[prepa[0]][2],
                                 'Enriched_Condition': enrichment_type_dict[prepa[1]],
                                 'Subproteome': subproteome_dict[prepa[1]]}
            else:
                preparation[prepa[1]]= {'Enrichment': enriched,
                                 'Dilution': prep[prepa[0]][0],
                                 'Added_vol': prep[prepa[0]][1],
                                 'Sample_vol': prep[prepa[0]][2],
                                 'Enriched_Condition': enrichment_type_dict[prepa[1]],
                                 'Subproteome': False}
        return preparation
    
    def parameter_gen(clean, params, conditions):

        N = len(conditions)
    
        template = pd.DataFrame(columns=params,
                     index=conditions)
    
        rng = np.random.default_rng(12345)
    
        for row in params:
    
            rand = rng.integers(low=0, high=10, size=N)
            if row == 'Cell count':
                rand = rand*1e8
            template[row] = rand
    
        template = template.replace(0, np.nan).reset_index().rename(columns={'index':'Conditions'})
        
        return template
   
    def generate_example_params(df):
        """
        Given a dataframe `df` containing columns 'Condition', 'Replicate', and 'Sample',
        generate an example parameter table with random values for the other columns.
        """
    
        # Get unique values for the 'Condition', 'Replicate', and 'Sample' columns
        conditions = df['Condition'].unique()
    
        # Generate random values for the other columns
        n_conditions = len(conditions)
    
        param_names = ['SampleVolume', 'ProteinConcentration', 'AmountMS',
                       'CellsPerML', 'TotalCultureVolume', 
                       'ProteinSRM', 'fmolSRM', 
                       'Enrichment', 'EnrichmentMode', 'StdDilution', 'StdVolume']
        param_types = [float, float, float,
                       float, float, 
                       str, float,
                       bool, str, float, float]
    
        # Generate random values for each parameter
        data = {}
        for name, dtype in zip(param_names, param_types):
            if dtype == bool:
                data[name] = np.random.choice([True, False], size=(n_conditions))
            elif name == 'EnrichmentMode':
                data[name] = np.random.choice(['Sampling', 'Amplification'], size=(n_conditions))
            elif name == 'ProteinSRM':
                data[name] = np.random.choice(df.Accession, size=(n_conditions))
            else:
                data[name] = np.random.rand(n_conditions) * 10
    
        # Create dataframe
        index = conditions
        df_params = pd.DataFrame(data=data, index=index, columns=param_names).reset_index(names='Condition')
    
        return df_params
 
    def where_to_find_std(clean, standards, thresh = 12):
    
        stan_list = standards['Accession'].unique()  
    
        where_to = clean[clean.Accession.isin(stan_list)].dropna(
                                ).groupby(['Sample'])['Accession'].nunique().reset_index()
    
        suggestion = where_to[where_to.Accession > thresh]['Sample'].to_list()
        
        return suggestion
        
    
    def pcomponent(clean, lfq_method='iBAQ', index=['Sample', 'Condition'], components=5, conditions=False):
        
        pivot_data = clean.dropna(subset=lfq_method).pivot_table(index=index, 
                                                         columns='Accession', 
                                                         values=lfq_method).fillna(0).reset_index()
        features = pivot_data.columns[2:]
        # Separating out the features
        x = pivot_data.loc[:, features].values
        # Separating out the target
        y = pivot_data.loc[:,['Sample']].values
        # Standardizing the features
        x = StandardScaler().fit_transform(x)
        pd.DataFrame(x, index=pivot_data.Sample)
        
        #pivot_data
        
        pca = PCA(n_components=components)
        principalComponents = pca.fit_transform(x)
       
        variance = pca.explained_variance_ratio_
        
        columns = [f'PC{var[0]+1}' for var in enumerate(variance)]
        
        principalDf = pd.DataFrame(data = principalComponents
                     , columns = columns)
        
        finalDF = principalDf.set_index(pivot_data.Sample).reset_index()
        finalDF['Condition'] = pivot_data.Condition
               
        
        return finalDF, columns, variance
    
    def create_random_df(df):
        # Select 10 random entries from the original data frame
        random_indices = np.random.choice(df.index, size=9, replace=False)
        selector = [col for col in df.columns if col in ["Accession", "Mol. weight [kDa]"]]
            
        random_df = df.loc[random_indices, selector]
        if "Mol. weight [kDa]" not in selector:
            random_df["Mol. weight [kDa]"] = np.random.rand(9)*100
    
        # Add fmol values of 50, 500, and 5000
        fmol_values = [50, 500, 5000] * 3
        random_df["StdConcentration"] = fmol_values
    
        return random_df
    
    def pivoter(df):
         
        cols = [col for col in df.select_dtypes(include=np.number).columns[1:]]
        values = [('mean', 'std') for i in range(len(cols))]
        
        operator = dict(zip(cols, values))
    
        df_mean = df.groupby(['Condition', 'Protein', 'Accession']).agg(operator,
        ).dropna().reset_index()
        
       # df_mean = df_mean[~df_mean.Accession.isin(standards)]
        
        df_mean.columns = [''.join(col).strip() for col in df_mean.columns.values]
        
        df_pivot = df_mean.pivot_table(index=['Protein', 'Accession'], columns='Condition',
                           values=df_mean.columns[2:]).swaplevel(axis=1)#.to_excel('Supp_table_X_secreted_proteins_diamide.xlsx')
        
        df_pivot.columns = [''.join(col).strip() for col in df_pivot.columns.values]
        
        cols = df_pivot.columns.to_list()#.sort()
        cols.sort()
        df_pivot = df_pivot[cols].reset_index()
        
        return df_pivot
        


