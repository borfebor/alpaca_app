import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer, SimpleImputer
import streamlit as st
from itertools import takewhile
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

class alpaca:
    
    def eats(file):
        try:
            df = pd.read_csv(file, sep='\t')
        except:
            try:
                df = pd.read_csv(file, sep=',')
            except:
                df = pd.read_excel(file)
            else:
                print('Sorry, this app ')

       	return df
       
    def eats_better(item):
        
        if 'txt' in item.name:
            
            df = pd.read_csv(item, sep='\t')
            
        elif 'csv' in item.name:

            df = pd.read_csv(item, sep=',')

        elif 'xlsx' in item.name:   
            
            df = pd.read_excel(item)
        
        else:
            
             st.warning('Warning: Could not read file.')

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
                        var_name='Sample', value_name=lfq_method)
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
        
    def spits(df, lfq_method='iBAQ', cleaning=True, formatting=True, identifier=None, normalization=False,
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
        
        default = ['Accession', 'Gene names', 'Mol. weight [kDa]']
        
        samples = [col for col in columns if lfq_method in col if '_' in col ]
        all_ids = [col for col in df.columns if col not in samples]
        
        wanted_ids = st.sidebar.multiselect('Data identifiers of interest', all_ids , default)
        ids = [col for col in df.columns if col in wanted_ids]
        conditions = list(set([item[len(lfq_method)+1:-3] for item in samples]))
        
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
            
            lfq = list(df.Sample.unique())
            prefix = ''.join(c[0] for c in takewhile(lambda x: all(x[0] == y for y in x), zip(*lfq)))
            df['Condition'] = df['Sample'].str[len(prefix):-3]
            df['Replicate'] = 'Replicate' + df['Sample'].str[-3:]
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


    def log_transform(df):
        ibaq = [x for x in df.columns if 'iBAQ' in x]
        df[ibaq] = np.log2(df[ibaq])
        return df

    def experimenter(df):
        """
        Looks for the experimental conditions in the samples and recognises how many (n)
        and which ones (condition)
        """
        try:
            condition = [x[5:-3] for x in df.columns if '_' in x]
            n = len(set(condition))  # defines the number of condition in the experiment
            r = int(len(condition) / n)
            condition = list(set(condition))
            print('Your experiment has', n, 'experimental conditions', condition, ', with', r, 'replicates.')
        except:
            print("I'm sorry, I could't predict your experimental conditions.")
        return n, r, condition

    def replicator(n):
        """
        Gets the amount of replicates from experimenter() and creates tags for the columns.
        It helps solving the problem of variable amount of replicates in different experiments.
        It returns a list of replicate names with numbers, that will be used later by condition_format() function.
        """
        i = 1
        rep = list()
        for x in range(n + 1):
            replicate = 'Replicate_' + str(i)
            rep.append(replicate)
            i += 1
        print('Generated your replicate names:', rep)
        return rep

    def condition_format(df, condition, rep):
        condition_col = list()
        condition_col = [x for x in df.columns if condition in x]
        try:
            cd = df[['Accession', 'Gene names', 'Mol. weight [kDa]']]
        except:
            print("Sorry, I couldn't found the columns Accession', 'Gene names', 'Mol. weight [kDa]")

        cd[rep] = df[condition_col]
        cd['Mean'] = cd[rep].median(axis=1)
        cd['Condition'] = condition
        cd = cd.dropna(subset=rep, thresh=2)
        return cd

    def sample_volumes(enrichment_type_dict, prep):
        volume_matcher = list()
        i = 2
        for key, value in enrichment_type_dict.items():
            volume_matcher.append([key, prep[i]])
            i += 3
        return volume_matcher
    
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
    
    def Quantify_old(df, condition, rep, coef, inter):
        appended = list()
        # This new loop solves the problem of variable amount of experimental conditions.
        # It creates a dataframe (clean_data) based on the condition list.
        for x in condition:
            cond = alpaca.condition_format(df, x, rep)
            appended.append(cond)

        clean_data = pd.concat(appended)
        clean_data['Amount (fmol)'] = 2 ** ((clean_data[['Mean']] - inter) / coef)
        return clean_data
    
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
    
    def wool(dummy, preparation, count_dict=None, enrichment_factors=None):  
        
        if enrichment_factors != None:
            dummy['Enriched_fmol'] = dummy.fmol.copy()
            for condition, enrichment in enrichment_factors.items():
                if len(condition) == 2:
                    print(f'{condition[1]} in subproteome {condition[0]} has been {round(enrichment, 2)} times enriched.')
                    dummy['Enriched_fmol'] = np.where((dummy.Condition == condition[1]) 
                                                      & (dummy.Subproteome == condition[0]),
                                                      dummy.fmol * enrichment, dummy.Enriched_fmol)
                elif len(condition) == 1:
                    print(f'{condition} has been {round(enrichment, 2)} times enriched.')
                    dummy['Enriched_fmol'] = np.where(dummy.Condition == condition, 
                                                      dummy.fmol * enrichment, dummy.Enriched_fmol)
                
            dummy['Molecules'] = dummy['Enriched_fmol'] * (
                            6.023e8)  # Avogadro's number fixed for fmol (-15)
        else:
            dummy['Molecules'] = dummy['fmol'] * (
                            6.023e8)  # Avogadro's number fixed for fmol (-15)
        
                    
        
        dummy['Molecules_per_cell'] = np.nan
        for prepo, values in preparation.items():
            if values['Subproteome'] == False:
                print('No subproteomic approach')
            for condition in values['Enriched_Condition']:
                
                vol = values['Sample_vol']
                cells = count_dict[condition] * vol
                
                if values['Subproteome'] != False:
                    for subprot in values['Subproteome']:
                        print(f'{condition} in {subprot} subproteome sample volume was {vol} ml.\nProteins quantified for {cells} cells.')
                        dummy['Molecules_per_cell'] = np.where((dummy.Condition == condition) & (dummy.Subproteome == subprot),
                                                               dummy.Molecules / cells, dummy.Molecules_per_cell)
                else:
                    print(f'{condition} sample volume was {vol} ml.\nProteins quantified for {cells} cells.')
                    dummy['Molecules_per_cell'] = np.where(dummy.Condition == condition,
                                                               dummy.Molecules / cells, dummy.Molecules_per_cell)
                    
        return dummy
    
    def wooler(df, preparation):

        enrichment_params = ['Enrichment', 'EnrichmentDirection', 'CorrectionFactorSRM', 'EnrichmentFactor']
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
                                        df.fmol / values['EnrichmentFactor'] * values['CorrectionFactorSRM'], 
                                        df.fmol)
                elif values['EnrichmentDirection'] == 'Down': 
                    """
                    This calculation is made for samples which correspond to a smaller fraction
                    to the original proteome. E.g., Secretome
                    """
                    df['fmol'] = np.where(df.Condition == condition, 
                                          df.fmol * values['EnrichmentFactor'] * values['CorrectionFactorSRM'], 
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
    
    def wool_old(dummy, prep, count_dict, conditions=None, enrichment_factors=None, enrichment_type_dict=None, subproteome=None):  
        
        if conditions == None:
            conditions = list(dummy.Condition.unique())
            
        if enrichment_type_dict == None:
            enrichment_type_dict = {'Enrichment_1' : list(dummy.Condition.unique())}
            
        if enrichment_factors != None:
            dummy['Enriched_fmol'] = dummy.fmol.copy()
            for condition, enrichment in enrichment_factors.items():
                if len(condition) == 2:
                    print(f'{condition[1]} in subproteome {condition[0]} has been {round(enrichment, 2)} times enriched.')
                    dummy['Enriched_fmol'] = np.where((dummy.Condition == condition[1]) & (dummy.Subproteome == condition[0]), dummy.fmol * enrichment, dummy.Enriched_fmol)
                elif len(condition) == 1:
                    print(f'{condition} has been {round(enrichment, 2)} times enriched.')
                    dummy['Enriched_fmol'] = np.where(dummy.Condition == condition, dummy.fmol * enrichment, dummy.Enriched_fmol)
                
            dummy['Molecules'] = dummy['Enriched_fmol'] * (
                            6.023e8)  # Avogadro's number fixed for fmol (-15)
        else:
            dummy['Molecules'] = dummy['fmol'] * (
                            6.023e8)  # Avogadro's number fixed for fmol (-15)
        
                    
        
        token = 0
        dummy['Molecules_per_cell'] = np.nan
        for item in enumerate(enrichment_type_dict.values()):
            for condition in enumerate(conditions):
                if condition[1] in item[1]:
                    accessor = 2
                    if token < item[0]:
                        accessor += 3
                    print(f'In {condition[1]}, there are {count_dict[condition[1]]} cells per ml in {prep[accessor]} ml')
                    cells = count_dict[condition[1]] * prep[accessor]
                    dummy['Molecules_per_cell'] = np.where((dummy.Condition == condition[1]) & (dummy.Subproteome == subproteome), dummy.Molecules / cells, dummy.Molecules_per_cell)
                    token = item[0]
            
        return dummy

    def std_amount(df, list, i):  # enrichments -> Integer. Enrichment types in our experiment, that's the reason that the function is in comprised in a loop
        e1 = df.copy()
        """
        This function will be comprised in the enricher() function, as a tool for calculating the enrichment amounts in
        every condition.
        """
        e1['Enrichment_type'] = 'Enrichment_' + str(i)
        e1['Added Standard proteins (ng)'] = e1[['Mix concentration (µg/µl)']]/ list[0] * list[1]
        e1['Amount (fmol) in sample'] = e1['Added Standard proteins (ng)'] / e1['MW (kDa)'] * 1000
        e1['Standard fmol/ml sample'] = e1['Amount (fmol) in sample'] / list[2]
        return e1

    def enricher(std, clean_data, prep, enrichment_type_dict, n, condition):
        # std -> dataframe including Enrichment standards
        # clean_data -> dataframe including our experimental data
        # prep - list with enrichment preparations (dilution, added volume, sample volume) * n
        # enrichments -> Integer. Enrichment types in our experiment
        # n -> Number of samples - output from experimenter
        """
        This function opens the possibility of introducing a variable number of enrichment conditions.
        It needs an extra step at the end for sorting the output, cause it calculates all the possible enrichments for
        all the standards. 
        """
        enrichments = len(enrichment_type_dict) 
        std_list = list()
        y = 0
        z = 3    
        x = 1
        for i in range(enrichments):
            tope = prep[y: z]
            cond = alpaca.std_amount(std, tope, x)
            std_list.append(cond)
            if x < enrichments+1:
                x +=1
            y += 3
            z += 3
        e1 = pd.concat(std_list)
        std_condition = list()
        for p in range(n):
            cond = clean_data[clean_data['Condition'] == condition[p]]
            std_condition.append(cond)
        std = pd.concat(std_condition)
        data = pd.merge(std, e1, on='Accession', how='inner')
        
        data = alpaca.enrichment_filter(data, enrichment_type_dict)
        with pd.ExcelWriter('enrichment_std.xlsx') as writer:  
            data.to_excel(writer, sheet_name='enrichment_std')
        return data

    def enrichment_filter(std_data, enrichment_type_dict):
        filtered = list()
        for key, value in enrichment_type_dict.items():
            for v in value:
                std_data_filtered = std_data[(std_data['Enrichment_type'] == key) & (std_data['Condition'] == v)]
                filtered.append(std_data_filtered)
        std_data = pd.concat(filtered)
        return std_data

    def enrichment_factors(data_std, enrichment_type_dict):
        enrichment_factors = dict(data_std.groupby(["Enrichment_type"])["Enrichment"].median())
        return enrichment_factors

    def getEnriched(dummy, enrichment_type_dict):
        enriqui_martin = list()
        for key, value in enrichment_type_dict.items():
            for v in value:
                dummy['Enrichment_type'] = key
                new_df = dummy[dummy.Condition == v]
                enriqui_martin.append(new_df)
        new_df = pd.concat(enriqui_martin, ignore_index=True)
        return new_df

    def getMolecules(dummy, enrichment_factors):
        enriqui_puig = list()

        for key, value in enrichment_factors.items():
            dummy['Enrichment_factor'] = value
            new_df = dummy[dummy.Enrichment_type == key]
            enriqui_puig.append(new_df)
        dummy = pd.concat(enriqui_puig)
        dummy['Enriched_fmol'] = dummy['Enrichment_factor'] * dummy['fmol']
        dummy['Molecules'] = dummy['Enriched_fmol'] * (
                    6.023e8)  # Avogadro's number fixed for fmol (-15)
        return dummy

    def getMolecules_cell(dummy, enrichment_type_dict, prep, count_dict):
        sample_volumes = alpaca.sample_volumes(enrichment_type_dict, prep)

        enriqui_viii = list()
        for k, v in count_dict.items():
            for key, value in enrichment_type_dict.items():
                for valor in value:
                    print(key, valor, v)
                    if k == valor:
                        dummy['Cells_per_ml'] = v
                        new_df = dummy.loc[(dummy['Condition'] == valor) & (dummy['Enrichment_type'] == key)]
                        enriqui_viii.append(new_df)
        la_tabla = pd.concat(enriqui_viii)

        dummy = list()
        for thing in sample_volumes:
            la_tabla['Molecules_per_cell'] = la_tabla['Molecules'] / (la_tabla['Cells_per_ml'] * thing[1])
            new_df = la_tabla.loc[la_tabla['Enrichment_type'] == thing[0]]
            dummy.append(new_df)
        la_tabla = pd.concat(dummy)
        return la_tabla
    
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

    @st.cache_data    
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
                       'CorrectionFactorSRM', 
                       'Enrichment', 'EnrichmentDirection', 'StdDilution', 'StdVolume']
        param_types = [float, float, float,
                       float, float, float,
                       bool, str, float, float]
    
        # Generate random values for each parameter
        data = {}
        for name, dtype in zip(param_names, param_types):
            if dtype == bool:
                data[name] = np.random.choice([True, False], size=(n_conditions))
            elif dtype == str:
                data[name] = np.random.choice(['Up', 'Down'], size=(n_conditions))
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
    
    def pregunton(df, query, look_on):
        
            if query.count(', ') >= 1:
                query = query.split(', ')
            elif query.count(',') >= 1:
                query = query.split(',')
            elif query.count(' ') >= 1:
                query = query.split(' ')
            else:
                query = [query]
            
                
            df = df[df[look_on].isin(query)]

            
            return df
            