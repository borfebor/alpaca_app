import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer, SimpleImputer
import streamlit as st
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from thefuzz import fuzz
import re
import os

class alpaca:
    
    def eats(file):
        
        if type(file) == str:
            file_name = file
        else:
            file_name = file.name
        
        if 'TXT' in file_name.upper():
            
            df = pd.read_csv(file, sep='\t')

        elif 'TSV' in file_name.upper():
            
            df = pd.read_csv(file, sep='\t')
            
        elif 'CSV' in file_name.upper():
            
            df = pd.read_csv(file, sep=',')
            
        elif 'XLSX' in file_name.upper():
                
            df = pd.read_excel(file)
        
        else:
            st.warning('Not compatible format')

       	return df
       
    def eats_better(file):
        
        file_type = file.name.split('.')[1].upper()
        
        if file_type == 'TXT':
            
            df = pd.read_csv(file, sep='\t')

        if file_type == 'TSV':
            
            df = pd.read_csv(file, sep='\t')
            
        if file_type == 'CSV':
            
            df = pd.read_csv(file, sep=',')
            
        if file_type == 'XLSX':
                
            df = pd.read_excel(file)
        
        else:
            st.warning('Not compatible format')

       	return df

    def data_cleaner(df, to_remove):
    
    	for col in to_remove:
        	if col in df.columns:
            		df = df.loc[lambda df: df[col].isna()]
    	df = df.drop(columns=to_remove)

    	return df

    
    def quant_norm(df):
        ranks = (df.rank(method="first")
                  .stack())
        rank_mean = (df.stack()
                       .groupby(ranks)
                       .mean())
        # Add interpolated values in between ranks
        finer_ranks = ((rank_mean.index+0.5).to_list() +
                        rank_mean.index.to_list())
        rank_mean = rank_mean.reindex(finer_ranks).sort_index().interpolate()
        return (df.rank(method='average')
                  .stack()
                  .map(rank_mean)
                  .unstack())
    
    
    def normalizer(data, lfq_method='iBAQ', normalization='Median', 
                   id_col=['Protein', 'Accession', 'Unique peptides', 'Mol. weight [kDa]']):
    
        df = data.copy()  
        print(f'Samples are normalized through {normalization} normalization')
        if normalization == 'Median':
    
            for sample in df.Sample.unique():
                operator = df[df.Sample == sample][lfq_method].median()
                
                df[lfq_method] = np.where(df.Sample == sample, df[lfq_method] - operator, df[lfq_method])
            new_col = f'm{lfq_method}'
            df = df.rename(columns={lfq_method:new_col})
                
        elif normalization == 'Relative':
            
            df[lfq_method] = np.power(df[lfq_method], 2)
            for sample in df.Sample.unique():
                
                operator = df[df.Sample == sample][lfq_method].sum()
                df[lfq_method] = np.where(df.Sample == sample, df[lfq_method] / operator, df[lfq_method])
            
            df[lfq_method] = np.log2(df[lfq_method])
                
            new_col = f'r{lfq_method}'
            df = df.rename(columns={lfq_method:new_col})
                
        elif normalization == 'Quantile':
    
            pivot_df = df.pivot_table(index=id_col, columns='Sample', values=lfq_method)
            
            lfq = [col for col in pivot_df.columns if lfq_method in col if '_' in col]
                
            pivot_df[lfq] = alpaca.quant_norm(pivot_df[lfq])
            
            df = pivot_df.reset_index().melt(id_vars=id_col, value_vars=lfq, var_name='Sample', value_name=lfq_method)
            new_col = f'q{lfq_method}'
            df = df.rename(columns={lfq_method:new_col})
            
        return df, new_col


    def machine_vision(df, conditions, ids, lfq_cols, lfq_method, identifier=None):
    	
        df_melt = df.melt(id_vars=ids, 
                        value_vars=lfq_cols,
                        var_name='Sample', value_name='value')
        df_melt = df_melt.rename(columns={'value':lfq_method})
        if identifier != None:
                df_melt = alpaca.identifiers(df_melt, identifier)
            
        return df_melt
    
    def identifiers(df, identifier):
        if type(identifier) is not dict:
            
            raise ValueError("A dictionary should be used to add identifiers (e.g. {'Subproteome':'Membrane'})")
        
        else:
            
            clean = df.copy()
            
            for col_name in identifier:
            
                if type(identifier[col_name]) is dict:
                    
                    for wish in identifier[col_name]:
            
                        column = [col for col in clean.columns if wish in clean[col].to_list()]
                        
                        #if col
                        if col_name not in clean.columns:
                            clean[col_name] = np.nan
                            
                        clean[col_name] = np.where(clean[column[0]] == wish,
                                                    identifier[col_name][wish], 
                                                    clean[col_name])
                elif type(identifier[col_name]) is str:
                    clean[col_name] = identifier[col_name]
        
            return clean
        
    def spits(df, lfq_method='iBAQ', cleaning=True, formatting=True, 
              lfq_columns=['iBAQ', 'LFQ', 'Top3', 'Intensity', 'MS/MS count'], 
              normalization=False, identifier=None,
              protein_ids=['Accession', 'Gene names', 'Mol. weight [kDa]']):
        '''
        

        Parameters
        ----------
        df : dataframe
            DESCRIPTION.
        lfq_method : str, optional
            DESCRIPTION. The default is 'iBAQ'.
        cleaning : bool, optional
            Removes potential contaminants, reverse and identified by site. The default is True.
        formating : bool, optional
            Rearranges the data with the desired structure for further analysis. The default is True.
        identifier : dict, optional
            Adds a column with the given identifiers (e.g. {'col_name':'value'}. The default is None.
        valid_values : int, optional
            Minimum quantified values per protein needed considering it as a valid quantification. The default is 2.
        protein_ids : list, optional
            List of headers which are desired to be used as identifiers. They should be present in the ProteinGroups.txt headers. 
            The default is ['Accession', 'Gene names', 'Unique peptides', 'Mol. weight [kDa]'].
        Returns
        -------
        df : TYPE
            DESCRIPTION.
        conditions : TYPE
            DESCRIPTION.

        '''
    
        df.columns = df.columns.str.replace('.: ', '')
        
        if 'Accession' not in df.columns:
            uniprot_key = [col for col in range(len(df.columns)) if 'Protein ID' in df.columns[col]]
            columns = list(df.columns)
            columns[uniprot_key[0]] = 'Accession'
            
            df.columns = columns
        
    	# Checking for data cleaning
        
        potential_cols = ['identified by site', 'contaminant', 'Reverse']
        cont_key = [col for col in df.columns for item in potential_cols if item in col]
        
        suggested = ['Accession', 'Gene names', 'Mol. weight [kDa]']
    
        conditions, samples, replicate_dict = alpaca.path_finder(df, lfq_method)
        all_ids = [col for col in df.columns if col not in samples]
        default = [col for col in all_ids if col in suggested]
        
        wanted_ids = st.sidebar.multiselect('Data identifiers of interest', all_ids , default)
        ids = [col for col in df.columns if col in wanted_ids]
        
        if cleaning is True:
            to_remove = st.sidebar.multiselect('Items to remove', cont_key, cont_key)
            df = alpaca.data_cleaner(df, to_remove)
            print(f'Items marked on {to_remove} have been removed from the dataset.')
            
        
        if formatting is True:
            
            df = alpaca.machine_vision(df, conditions, ids, samples, lfq_method, identifier)
            
            df = df.replace(0, np.nan)
            df[lfq_method] = np.log2(df[lfq_method])
            
            if normalization != False:
                
                df, lfq_method = alpaca.normalizer(df, lfq_method=lfq_method, 
                                                   normalization=normalization, id_col=ids)
                
            df['Condition'] = np.nan
            df['Replicate'] = np.nan
            
            for item in replicate_dict:
                
                df['Condition'] = np.where(df.Sample == item, replicate_dict[item][0], df.Condition)
                df['Replicate'] = np.where(df.Sample == item, replicate_dict[item][1], df.Replicate)

            df['Sample'] = df['Sample'].str.rsplit(' ', expand=True, n=1)[1]
            df['Condition'] = df['Condition'].str.replace(r'[_-]0', '', regex=True)

            df = df.dropna(subset=lfq_method)
            
            if 'Gene names' in ids:
                df['Gene names'] = df['Gene names'].str[0].str.upper() + df['Gene names'].str[1:]
                df = df.rename(columns={'Gene names':'Protein'})
            print('Dataset formated for further analysis and visualisation.')
            
            conditions = df.Condition.unique()
        
        else:
            ids = ids + samples
            df = df[ids]
            print('Data is formated for human vision.\nThat could lead to errors or incompatibilities in further analysis using Alpaca pipeline.\nConsider formating your dataset if you see any anomally.')
            
        return df, conditions, lfq_method    
    
    def scientist(prep, enrichment_type_dict, subproteome_dict=None):
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
    
    def abacus(ups2, concentration=0.5, in_sample=6.0, total_protein=10):
        
        #ups2 = alpaca.eats(standards)
        
        fmol_col = [fmol for fmol in ups2.columns if ('mol' or 'MOL') in fmol]
        MW = [mw for mw in ups2.columns if ('Da' or 'MW') in mw]
        
        print('Got column:', fmol_col[0], '- Calculating fmols for you')
        µg_standards = in_sample * concentration
        
        ups2['fmol_inSample'] = ups2[fmol_col[0]] / 10.6 * µg_standards
        ups2['log2_Amount_fmol'] = np.log2(ups2.fmol_inSample)
        ups2['Mass fraction (fmol/µg_extract)'] = ups2.fmol_inSample * (µg_standards / total_protein)
        ups2['log2_Mass_fract'] = np.log2(ups2['Mass fraction (fmol/µg_extract)'])
        
        volume = 10.6 / concentration  # concentration in µL
        
        print('UPS2 standards vial concentration:', concentration, 'µg/µl | Resuspended in:', volume, 'µl')
        ups2['Stock_ng_µl'] = ups2[fmol_col[0]] * ups2[MW[0]] / (volume*1e6)
        
        added = in_sample * ups2['Stock_ng_µl']
        print(in_sample, 'µl added to the sample')
        ups2['In_sample_ng'] = ups2['Stock_ng_µl'] * added
        return ups2
    
    def regression(df, ups2, lfq_col='iBAQ', filter_col='Replicate', added_samples=None, valid_values=2):
    
        data = pd.merge(ups2, df, on='Accession', how='right')
        #data = data.dropna(subset=[fmol_col])
        if added_samples != None:
            data = data[data[filter_col].isin(added_samples)]
            
        ups_red = data.dropna(subset=lfq_col).groupby(['Accession', 'log2_Amount_fmol']).apply(lambda x: pd.Series({
                                lfq_col: x[lfq_col].median(), 'N': x[lfq_col].nunique()})).reset_index()
        
        ups_red = ups_red[ups_red.N >= valid_values] #Filters for a minimum of values needed for fitting
        
        X = ups_red['log2_Amount_fmol'].values.reshape(-1, 1)  # values converts it into a numpy array
        Y = ups_red[lfq_col].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
        linear_regressor = LinearRegression().fit(X, Y)  # create object for the class & perform linear regression
        Y_pred = linear_regressor.predict(X)  # make predictions
        # The coefficients
        coef = linear_regressor.coef_
        inter = linear_regressor.intercept_
        print(f'Coefficients: {coef}')
        print(f'Intercept: {inter}')
        # The mean squared error
        print('Mean squared error: %.2f'
            % mean_squared_error(Y, Y_pred))
        # The coefficient of determination: 1 is perfect prediction
        R2 = r2_score(Y, Y_pred)
        print('Coefficient of determination: %.2f'
            % R2)
        
        return ups_red, coef[0], inter, R2
    
    def abacusreg(ups_red, lfq_col='iBAQ', R2='', save=True):
        
        sns.set_context('paper')
        sns.set(style='whitegrid', font_scale=1.5)
        g = sns.lmplot(x='log2_Amount_fmol', y=lfq_col, data=ups_red, palette=['#656565'],  aspect=1)
        g.set_axis_labels("UPS2 amount (log2)", f"UPS2 {lfq_col} (log2)")
        plt.text(ups_red['log2_Amount_fmol'].min()-0.2, ups_red[lfq_col].max()-0.5, round(R2, 3), fontsize=20)
        if save == True:
            g.savefig('UPS2_Quantification.svg', bbox_inches='tight', pad_inches=0.5)
            
    def moles(df, coef, inter, lfq_col='iBAQ', ratio= 1, ratio_col='identifier'):
    
        df['fmol'] = 2 ** ((df[lfq_col] - inter) / coef)
        
        if type(ratio) is int:
            df['fmol'] = df['fmol'] / ratio
        elif type(ratio) is dict: 
            for key, value in ratio.items():
                df['fmol'] = np.where(df[ratio_col] == key,  df['fmol'] / value,  df['fmol'])
        
        return df
        
    def census(df, standards, concentration=0.5, in_sample=6.0, lfq_col='iBAQ', ratio=1, 
               total_protein= 1, filter_col='Replicate', added_samples=None, valid_values=2, save=True):
        '''
        

        Parameters
        ----------
        df : Dataframe
            Clean data from quantified proteins
        standards : File (.csv, .txt, .xlsx)
            UPS2 dynamic standards information.
        concentration : float, optional
            Standards stock concentration. The default is 0.5.
        in_sample : float, optional
            Added volume in sample. The default is 6.0.
        lfq_col : ('iBAQ', 'LFQ', 'Intensity'), optional
            DESCRIPTION. The default is 'iBAQ'.
        added_samples : str or None, optional
            Samples (conditions) in which it was added the standards. The default is None.
        save : bool, optional
            Toggle save the graphs. The default is True.

        Returns
        -------
        df : dataframe
            Quantified proteins.
        ups_red : dataframe
            Measured standards in the sample.
        coef : float
            Regression slope.
        inter : float
            Regression interception.

        '''
        ups2 = alpaca.abacus(standards, concentration, in_sample, total_protein=total_protein) # Arranges the standards
        ups_red, coef, inter, R2 = alpaca.regression(df, ups2, lfq_col=lfq_col, filter_col=filter_col,
                                                     added_samples=added_samples, valid_values=valid_values) # Regression between intensities and standards
        alpaca.abacusreg(ups_red, lfq_col=lfq_col, R2=R2, save=save) # Plots the Regression
        df = alpaca.moles(df, coef, inter, lfq_col=lfq_col, ratio = ratio)
    
        return df, ups_red, coef, inter, R2
    
    def gathers(df, enrichment_standards, preparation, subproteome=None, QC=False, deviation_lim=10, thresh=10, 
                imputation=True, strategy='KNN', plot=False, save_plot=False, lfq_method='iBAQ'):
        '''
        

        Parameters
        ----------
        df : dataframe
            Quantified data.
        enrichment_standards : dataframe
            Enrichment standards data.
        preparation : dict, optional
            Dictionary with the gathering which samples are prepared in which way.
            If None, it will assume that all conditions were prepared the same way.
            The default is None.
        subproteome: str or list. The default is None.
        Subproteomic fractions that have been enriched
        QC : bool, optional
            Toggle for outlier QC. The default is True.
        deviation_lim : int or float, optional
            Parameter for detecting outliers. The default is 10.
        thresh: float or int, optional
            Parameter for filtering outlier values regarding the deviation to the enrichment median of a protein. The default is 10.
        imputation : bool, optional
            Toggle for imputation of dropped outliers. The default is True.
        strategy : str, optional
            Imputation method ('KNN', 'mean', 'median', 'most frequent', 'constant').. The default is 'KNN'.

        Returns
        -------
        e_test : df
            Dataframe containing the quantified spiked_in standards.
        enrichment_factors : dict
            Dict with our calculated enrichment factors.

        '''
            
        e_test = alpaca.enrichment_calculator(df, enrichment_standards, preparation, lfq_method).dropna(subset='Enrichment')
            
        if QC == True:
            e_test = alpaca.looksafter(e_test, deviation_lim, thresh, imputation, strategy)
            
        grouping = ['Condition', 'Replicate']
        col_grouper = [columns for columns in df.columns if columns in grouping]

        enrichments = e_test.groupby(col_grouper)['Enrichment'].median().reset_index()
        col_grouper.remove('Replicate')
        enrichment_factors = enrichments.groupby(col_grouper).apply(lambda x: pd.Series({
                                                                    'EnrichmentFactor':x["Enrichment"].median()})
            ).reset_index()

        preparation = preparation.merge(enrichment_factors, on='Condition', how='left')
            
        for key, value in enrichment_factors.iterrows():
            print(f'Enrichment factor on condition: {key} = {value}')
            
        if plot == True:
            sns.catplot(data=enrichments, x='Condition', y='Enrichment', kind='box', width=0.5)
        if save_plot == True:
            plt.savefig('spiked_in_standards.svg', bbox_inches='tight', pad_inches=0.5)
        
        return e_test, preparation
   
    def correctionSRM(df, preparation):

        data = pd.DataFrame()

        for index, content in preparation.iterrows():
            
                vals = str(content['fmolSRM']).split(',')
                values = [float(val) for val in vals]
    
                entries = content['ProteinSRM'].split(',')
    
                condition = content.Condition
    
                temp = pd.DataFrame(zip(entries,values), columns=['Accession', 'fmolSRM'])
                temp['Condition'] = condition
                data = pd.concat([data, temp])
                
        if data.shape[0] >= 1:
    
            correction = data.merge(df, 
                                   on=['Accession', 'Condition'],
                                   how='left').groupby([
                            'Condition', 'Accession']).apply(lambda x: pd.Series({
                            'CorrectionSRM': x['fmolSRM'].mean() / x['fmol'].mean()
                        })).reset_index(
            ).groupby('Condition')['CorrectionSRM'].mean().reset_index()

            preparation = preparation.merge(correction, on='Condition') 
    
        return preparation
        
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
    
    def standards_preparation_check(preparation):
    
        try:
            vol_col = [col for col in preparation.columns if ("STDV" or "STD_V" or "STD V") in col.upper()][0]
        except:
            print('Added volume is missing in sample parameters')
        
        try:
            splv_col = [col for col in preparation.columns if ("SAMPLEV" or "SAMPLE_V" or "SAMPLE V") in col.upper()][0]
        except:
            print('Sample volume is missing in sample parameters')
        
        try:
            dil_col = [col for col in preparation.columns if ("STDD" or "STD_D" or "STD D") in col.upper()][0]
        except:
            print('Dilution of the standard stock solution is missing in sample parameters')
            
        return vol_col, splv_col, dil_col
    
    def multiplier(enriched_conditions, standards):
    
        arranger = pd.DataFrame()
    
        for condition in enriched_conditions:
    
            standards['Condition'] = condition
            arranger = pd.concat([arranger, standards])
            
        return arranger
    
    def enrichment_calculator(df, standards, preparation, lfq_method='iBAQ', subproteome=None):
        
        enriched_conditions = [item[1]['Condition'] for item in preparation.iterrows() if item[1]['Enrichment'] == True]
        
        if enriched_conditions == []:
            pass
        else:
            vol_col, splv_col, dil_col = alpaca.standards_preparation_check(preparation)
    
            standards_mod = alpaca.multiplier(enriched_conditions, standards)
    
            standards_mod = standards_mod.merge(preparation, on=['Condition'])
            
            MW = [col for col in standards_mod.columns if 'Da' in col][0]
    
            standards_mod['StdMass'] = standards_mod['StdConcentration'] / standards_mod[dil_col] * standards_mod[vol_col]
            standards_mod['StdFmol'] = standards_mod['StdMass'] / standards_mod[MW] * 1000
    
            standards_mod['StdFmol'] = np.where(standards_mod['StdFmol'] == 0, np.nan, standards_mod['StdFmol'])
    
            standards_mod = standards_mod[['Accession', 'Condition', MW, 'StdMass','StdFmol']]
    
            ID_standards = df[df.Condition.isin(enriched_conditions)
                             ].merge(standards_mod, how='right', on=["Accession", "Condition"])
            
            ID_standards['Enrichment'] = ID_standards['StdFmol'] / ID_standards['fmol']  
    
            e_reduction = ['Accession', MW, 'Sample', lfq_method, 
                                           'Condition', 'Replicate', 'fmol', 'StdFmol', 'Subproteome', 'Enrichment']
    
            ID_standards = ID_standards[[col for col in ID_standards.columns if col in e_reduction]]
                    
            return ID_standards

    
    def droper(e_test, std=10, thresh=10):

        test = e_test.copy()
        test = test.dropna().reset_index(drop=True)
        for item in test.Accession.unique():
            for condition in test.Condition.unique():
                ef = test[(test.Accession == item) & (test.Condition == condition)]['Enrichment'].dropna()
                
                if ef.std() > std:
                    dropers = [item for item in np.array(ef) if np.absolute((item - ef.median())) > thresh]
                    for droper in dropers:
                        test['Enrichment'] = np.where(test.Enrichment == droper, np.nan, test.Enrichment)
                        new_enrichment = test[(test.Accession == item) &(test.Condition == condition)]['Enrichment']
                        new_std = round(new_enrichment.std(),2)
                    if len(new_enrichment.dropna()) < 2:
                        test = test.drop(index=ef.index)
                        print('Droped', item, 'at', condition, 'cause there was only 1 value.')
                    else:
                        print(f"Droped a replicate of {item} at {condition} cause STD was too high ({round(ef.std(),2)}).\nStandard deviation is now to {new_std}.")
                
                elif len(ef) < 2:
                    test = test.drop(index=ef.index)
                    print('Droped', item, 'at', condition, 'cause there was only 1 value.')
        return test

    def KNNimputation(test):
        test = test.reset_index(drop=True)
        imputed = list()
        for item in test.Accession.unique():
            for condition in test.Condition.unique():
                sorted_data = test[(test.Accession == item) & (test.Condition == condition)]
                ef = np.array(sorted_data.Enrichment)
                ef = ef.reshape(-1, 1) 
                if len(ef) >= 1:
                    imputer = KNNImputer(n_neighbors=len(ef))
                    ef_mod = imputer.fit_transform(ef)
                    ef_mod_std = ef.std()
                    ef_mod = np.reshape(ef_mod, len(ef_mod))
                    data = pd.DataFrame(ef_mod, index=sorted_data.index)
                    imputed.append(data)
            print(f'Enrichment outliers imputed through k-Nearest Neighbors.\nFor {item} at {condition}, standard deviation is now {round(ef_mod_std,2)}')
        test['Enrichment'] = pd.concat(imputed)
        return test
    
    def SimpleImputation(test, strategy='Median'):
        test = test.reset_index(drop=True)
        imputed = list()
        for item in test.Accession.unique():
            for condition in test.Condition.unique():
                sorted_data = test[(test.Accession == item) & (test.Condition == condition)]
                ef = np.array(sorted_data.Enrichment)
                ef = ef.reshape(-1, 1) 
                if len(ef) >= 1:
                    imputer = SimpleImputer(missing_values=np.nan, strategy=strategy)
                    ef_mod = imputer.fit_transform(ef)
                    ef_mod_std = ef.std()
                    ef_mod = np.reshape(ef_mod, len(ef_mod))
                    data = pd.DataFrame(ef_mod, index=sorted_data.index)
                    imputed.append(data)
        
            print(f'Enrichment outliers imputed through {strategy}.\nFor {item} at {condition}, standard deviation is now {round(ef_mod_std,2)}')
        test['Enrichment'] = pd.concat(imputed)
        return test
    
    def looksafter(test, deviation_lim=10, thresh=10, imputation=True, strategy='KNN'):
        '''
        

        Parameters
        ----------
        test : dataFrame
            Enrichment standards data.
        deviation_lim : int or float, optional
            Standard deviation threshold among a standard to drop a replicate. The default is 10.
        imputation : Bool, optional
            Option to enable imputation of values with high deviation. If False, outliers will be dropped. The default is True.
        strategy : str, optional
            Imputation method ('KNN', 'mean', 'median', 'most frequent', 'constant'). The default is 'KNN'.

        Returns
        -------
        test_imputed : dataFrame
            Clean enrichment standards data.

        '''
        test_clean = alpaca.droper(test, std=deviation_lim, thresh=thresh)
        if imputation == True:
            if strategy != 'KNN':
                test_imputed = alpaca.SimpleImputation(test_clean, strategy=strategy)
            else:
                test_imputed = alpaca.KNNimputation(test_clean)
        else:
            test_imputed = test_clean.copy()
        return test_imputed
    
    def wooler(df, preparation):

        enrichment_params = ['Enrichment', 'EnrichmentDirection', 'ProteinSRM', 'fmolSRM', 'EnrichmentFactor']
        sample_params = ['SampleVolume', 'ProteinConcentration', 'AmountMS']
        cells_params = ['CellsPerML', 'TotalCultureVolume']
        
        if 'EnrichmentFactor' in preparation.columns.to_list():
        
            for condition, values in preparation.set_index('Condition')[enrichment_params].fillna(1).iterrows():
                
                if values['EnrichmentDirection'] == 'Up':
                    """
                    This calculation is made for samples which correspond to a higher fraction 
                    compared to the original proteome. E.g., Membrane
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                        df.fmol / values['EnrichmentFactor'], 
                                        df.fmol)
                    
                elif values['EnrichmentDirection'] == 'Down': 
                    """
                    This calculation is made for samples which correspond to a smaller fraction
                    to the original proteome. E.g., Secretome
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                          df.fmol * values['EnrichmentFactor'], 
                                          df.fmol)
                    
            preparation = alpaca.correctionSRM(df, preparation)

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
                       'Enrichment', 'EnrichmentDirection', 'StdDilution', 'StdVolume']
        param_types = [float, float, float,
                       float, float, 
                       str, float,
                       bool, str, float, float]
    
        # Generate random values for each parameter
        data = {}
        for name, dtype in zip(param_names, param_types):
            if dtype == bool:
                data[name] = np.random.choice([True, False], size=(n_conditions))
            elif name == 'EnrichmentDirection':
                data[name] = np.random.choice(['Up', 'Down'], size=(n_conditions))
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
        
    def pca(clean, standards, lfq_method='iBAQ'):
        
        stan_list = standards['Accession'].unique()  
    
        idstd = clean[clean.Accession.isin(stan_list)].dropna()
        
        pivot_data = idstd.dropna(subset=lfq_method).pivot_table(index='Sample', 
                                                         columns='Accession', 
                                                         values=lfq_method).fillna(0).reset_index()
        features = pivot_data.columns[1:]
        # Separating out the features
        x = pivot_data.loc[:, features].values
        # Separating out the target
        y = pivot_data.loc[:,['Sample']].values
        # Standardizing the features
        x = StandardScaler().fit_transform(x)
        pd.DataFrame(x, index=pivot_data.Sample)
        #pivot_data
        
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(x)
        principalDf = pd.DataFrame(data = principalComponents
                     , columns = ['PC1', 'PC2'])
        
        finalDF = principalDf.set_index(pivot_data.Sample).sort_values(by='PC1', ascending=False).reset_index()
        
        return finalDF
    
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
    
    def KMeans(sorted_data, component='PC1', n_clusters=2):
    
        model = KMeans(n_clusters=n_clusters)
    
        model.fit(sorted_data[component].values.reshape([-1, 1]) )
    
        all_predictions = model.predict(sorted_data[component].values.reshape([-1, 1]))
        sorted_data['cluster'] = pd.Series(all_predictions)
        
        group = all_predictions[0]
        suggestions = sorted_data[sorted_data.cluster == group]['Sample']
        
        return sorted_data, suggestions
    
    def create_random_df(df):
        # Select 10 random entries from the original data frame
        random_indices = np.random.choice(df.index, size=9, replace=False)
        random_df = df.loc[random_indices, ["Accession", "Mol. weight [kDa]"]]
    
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
        
    def match_names(name, df, thresh=75):
    
        for index, x in enumerate(df.columns):

            score = fuzz.ratio(name.lower(), x.lower())

            if score > thresh:

                return index, name, score
        
    def matchmaker(df_original, df_desired):
    
        editable = df_original.columns.to_list()

        for x in df_desired.columns:

            match = alpaca.match_names(x, df_original, 75)

            editable[match[0]] = match[1]

        df_original.columns = editable

        return df_original
    
    def path_finder(df, lfq_method):

        candidates = [re.findall(r"\w+", col)[-1] for col in df.columns if lfq_method in col if len(col) > len(lfq_method)]

        conditions = []
        replicate = dict()
        previous = ''

        for candidate in permutations(candidates, 2):

            name1 = re.sub('[^0-9a-zA-Z]+', '', candidate[0])

            name2 = re.sub('[^0-9a-zA-Z]+', '', candidate[1])

            ratio = fuzz.ratio(name1, name2)

            if ratio > 85:

                common = os.path.commonprefix(candidate)

                if common in previous:

                    common = previous

                previous = common

                if common != '':

                    conditions.append(common)

                    rep = candidate[0].replace(common, '')

                    col = [col for col in df.columns for item in conditions if lfq_method in col if candidate[0] in col][0]

                    replicate[col] = [common, f'Replicate_{rep}']

        conditions = list(dict.fromkeys(conditions))
        sample_cols = list(replicate.keys())

        return conditions, sample_cols,  replicate

