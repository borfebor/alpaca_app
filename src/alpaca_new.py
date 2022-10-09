import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.impute import KNNImputer, SimpleImputer
import streamlit as st


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

    def data_cleaner(df, to_remove):
    
    	for col in to_remove:
        	if col in df.columns:
            		df = df.loc[lambda df: df[col].isna()]
    	df = df.drop(columns=to_remove)

    	return df


    def machine_vision(df, conditions, ids, lfq_cols, lfq_method, subproteome=None):
    	
        df_melt = df.melt(id_vars=ids, 
                        value_vars=lfq_cols,
                        var_name='Sample', value_name=lfq_method)
        if subproteome != None:
                df_melt['Subproteome'] = subproteome
        return df_melt
    
    def spits(df, lfq_method='iBAQ', cleaning=True, formating=True, subproteome=None):
        '''
        

        Parameters
        ----------
        df : TYPE
            DESCRIPTION.
        lfq_method : TYPE, optional
            DESCRIPTION. The default is 'iBAQ'.
        cleaning : TYPE, optional
            DESCRIPTION. The default is True.
        formating : TYPE, optional
            DESCRIPTION. The default is True.
        subproteome : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        df : TYPE
            DESCRIPTION.
        conditions : TYPE
            DESCRIPTION.

        '''
    
        df.columns = df.columns.str.replace('.: ', '')

        uniprot_key = [col for col in range(len(df.columns)) if 'Protein ID' in df.columns[col]]
        columns = list(df.columns)
        columns[uniprot_key[0]] = 'Accession'
            
        df.columns = columns

    	# Checking for data cleaning
        potential_cols = ['identified by site', 'contaminant', 'Reverse']
        cont_key = [col for col in df.columns for item in potential_cols if item in col]
        st.write(cont_key)
        
        default = ['Accession', 'Gene names', 'Unique peptides', 'Mol. weight [kDa]']
        
        samples = [col for col in columns if lfq_method in col if '_' in col ]
        
        wanted_ids = st.sidebar.multiselect('Data identifiers of interest', columns.remove(samples) , default)
        ids = [col for col in df.columns if col in wanted_ids]
        conditions = list(set([item[len(lfq_method)+1:-3] for item in samples]))
        
        if cleaning == True:
            to_remove = st.sidebar.multiselect('Items to remove', cont_key, cont_key)
            df = alpaca.data_cleaner(df, to_remove)
            print(f'Items marked on {to_remove} have been removed from the dataset.')
        
        if formating == True:
            df = alpaca.machine_vision(df, conditions, ids, samples, lfq_method, subproteome)
            df['Condition'] = df['Sample'].str[len(lfq_method)+1:-3]
            df['Replicate'] = 'Replicate' + df['Sample'].str[-3:]
            df = df.replace(0, np.nan)
            df[lfq_method] = np.log2(df[lfq_method])
            if 'Gene names' in ids:
                df['Gene names'] = df['Gene names'].str[0].str.upper() + df['Gene names'].str[1:]
                df = df.rename(columns={'Gene names':'Protein'})
            print('Dataset formated for further analysis and visualisation.')
        else:
            print('Data is formated for human vision.\nThat could lead to errors or incompatibilities in further analysis using Alpaca pipeline.\nConsider formating your dataset if you see any anomally.')
   
        return df, conditions
    
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
    
    def abacus(standards, concentration=0.5, in_sample=6.0):
        
        ups2 = alpaca.eats(standards)
        
        fmol_col = [fmol for fmol in ups2.columns if ('mol' or 'MOL') in fmol]
        MW = [mw for mw in ups2.columns if ('Da' or 'MW') in mw]
        
        print('Got column:',fmol_col[0], '- Calculating fmols for you')
        ups2['log2_Amount_fmol'] = np.log2(ups2[fmol_col[0]])
        ups2['Mass fraction (fmol/µg_extract)'] = ups2.log2_Amount_fmol / 10.6
        ups2['log2_Mass_fract'] = np.log2(ups2['Mass fraction (fmol/µg_extract)'])
        
        volume = 10.6 / concentration  # concentration in µL
        
        print('UPS2 standards vial concentration:', concentration, 'µg/µl | Resuspended in:', volume, 'µl')
        ups2['Stock_ng_µl'] = ups2[fmol_col[0]] * ups2[MW[0]] / (volume*1e6)
        
        added = in_sample * ups2['Stock_ng_µl']
        print(in_sample, 'µl added to the sample')
        ups2['In_sample_ng'] = ups2['Stock_ng_µl'] * added
        return ups2
    
    def regression(df, ups2, lfq_col='iBAQ', replicate=None):
    
        data = pd.merge(ups2, df, on='Accession', how='right')
        #data = data.dropna(subset=[fmol_col])
        if replicate != None:
            data = data[data['Replicate'] == replicate]
        ups_red = data.groupby(['Accession', 'log2_Amount_fmol'])[lfq_col].median().reset_index().dropna()
        
        X = ups_red['log2_Amount_fmol'].values.reshape(-1, 1)  # values converts it into a numpy array
        Y = ups_red[lfq_col].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
        linear_regressor = LinearRegression().fit(X, Y)  # create object for the class & perform linear regression
        Y_pred = linear_regressor.predict(X)  # make predictions
        # The coefficients
        coef = linear_regressor.coef_
        inter = linear_regressor.intercept_
        print('Coefficients: \n', coef)
        print('Intercept: \n', inter)
        # The mean squared error
        print('Mean squared error: %.2f'
            % mean_squared_error(Y, Y_pred))
        # The coefficient of determination: 1 is perfect prediction
        R2 = r2_score(Y, Y_pred)
        print('Coefficient of determination: %.2f'
            % R2)
        
        return ups_red, coef, inter, R2
    
    def abacusreg(ups_red, lfq_col='iBAQ', R2='', save=True):
        
        sns.set_context('paper')
        sns.set(style='whitegrid', font_scale=1.5)
        g = sns.lmplot(x='log2_Amount_fmol', y=lfq_col, data=ups_red, palette=['#656565'],  aspect=1)
        g.set_axis_labels("UPS2 amount (log2)", "UPS2 IBAQ (log2)")
        plt.text(ups_red['log2_Amount_fmol'].min()-0.2, ups_red[lfq_col].max()-0.5, round(R2, 3), fontsize=20)
        if save == True:
            g.savefig('UPS2_Quantification.svg', bbox_inches='tight', pad_inches=0.5)
            
    def moles(df, coef, inter, lfq_col='iBAQ',):
    
        df['fmol'] = 2 ** ((df[['iBAQ']] - inter) / coef)
        
        return df
        
    def census(df, standards, concentration=0.5, in_sample=6.0, lfq_col='iBAQ', replicate=None, save=True):
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
        replicate : str or None, optional
            Replicate in which it was added the standards. The default is None.
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
        ups2 = alpaca.abacus(standards, concentration, in_sample) # Arranges the standards
        ups_red, coef, inter, R2 = alpaca.regression(df, ups2, lfq_col=lfq_col, replicate=replicate) # Regression between intensities and standards
        alpaca.abacusreg(ups_red, lfq_col=lfq_col, R2=R2, save=save) # Plots the Regression
        df = alpaca.moles(df, coef, inter, lfq_col=lfq_col)
    
        return df, ups_red, coef, inter
    
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
    
    def gathers(df, enrichment_standards, preparation, subproteome=None, QC=True, deviation_lim=10, thresh=10, 
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
            

            
        e_test = alpaca.spiked_in(df, enrichment_standards, preparation, lfq_method, subproteome)
            
        if QC == True:
            e_test = alpaca.looksafter(e_test, deviation_lim, thresh, imputation, strategy)
            
        grouping = ['Condition', 'Subproteome', 'Replicate']
        col_grouper = [columns for columns in df.columns if columns in grouping]

        enrichments = e_test.groupby(col_grouper)['Enrichment'].median().reset_index()
        col_grouper.remove('Replicate')
        enrichment_factors = dict(enrichments.groupby(col_grouper)["Enrichment"].median())
            
        for key, value in enrichment_factors.items():
            print(f'Enrichment factor on condition: {key} = {value}')
            
        if plot == True:
            sns.catplot(data=enrichments, x='Condition', y='Enrichment', kind='box', width=0.5)
        if save_plot == True:
            plt.savefig('spiked_in_standards.svg', bbox_inches='tight', pad_inches=0.5)
        
        
        return e_test, enrichment_factors
    
        

    def spiked_in(clean_data, enrich, preparation, lfq_method='iBAQ', subproteome=None):
        
        enriched_conditions = [item['Enriched_Condition'] for item in preparation.values()]
        enrichment_dict = dict(zip(preparation.keys(), np.array_split(np.concatenate(enriched_conditions), len(preparation))))
        
        volumes = [[item['Dilution'], item['Added_vol'], item['Sample_vol']] for item in preparation.values()]
        volume_dict = dict(zip(preparation.keys(), np.array_split(np.concatenate(volumes), len(preparation))))
        lab_work = ['dilution', 'added', 'volume']
        
        e_match = pd.DataFrame(enrichment_dict).melt().drop_duplicates().rename(
                                    columns={'variable':'e_type', 'value':'Condition'})
        e_match['preparation'] = e_match.e_type.map(volume_dict)
        e_match[lab_work] = e_match.preparation.apply(pd.Series)
            
        e_test = clean_data.merge(enrich, on='Accession', 
                                       how='right').merge(e_match, on='Condition', how='right')
        if subproteome!=None:
                            
            if type(subproteome) == str:
                my_item = subproteome
                subproteome = list()
                subproteome.append(my_item)
                        
           
            e_test = e_test[e_test.Subproteome.isin(subproteome)] 
                
        e_test['Added_ng'] = e_test['Mix concentration (µg/µl)'] / e_test['dilution'] * e_test['added']
        e_test['Added_fmol'] = e_test['Added_ng'] / e_test['Mol. weight [kDa]'] * 1000
        
        e_test['Added_fmol'] = np.where(e_test['Added_fmol'] == 0, np.nan, e_test['Added_fmol'])
                
        e_reduction = ['Accession', 'Mol. weight [kDa]', 'Sample', lfq_method, 
                               'Condition', 'Replicate', 'fmol','Added_fmol', 'Subproteome']
                
        e_test = e_test[[col for col in e_test.columns if col in e_reduction]]
                
        e_test['Enrichment'] = e_test['fmol'] / e_test['Added_fmol'] 
                    
        return e_test
    
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
    
    def wool(dummy, preparation, count_dict, enrichment_factors=None):  
        
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
