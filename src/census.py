#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 18:34:39 2024

@author: borfebor
"""

import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import re
from src.spits import Normalization, tools, clean, Imputation, detective


class Quantification:

    def abacus(ups2, concentration=0.5, in_sample=6.0, total_protein=10):
        
        #ups2 = alpaca.eats(standards)
        
        fmol_col = [fmol for fmol in ups2.columns if ('mol' or 'MOL') in fmol]
       # MW = [mw for mw in ups2.columns if ('Da' or 'MW') in mw]
        
        print('Got column:', fmol_col[0], '- Calculating fmols for you')
        µg_standards = in_sample * concentration
        
        ups2['fmol_inSample'] = ups2[fmol_col[0]] / 10.6 * µg_standards
        ups2['log2_Amount_fmol'] = np.log2(ups2.fmol_inSample)

        return ups2
    
    def regression(df, standards, lfq_col='iBAQ', filter_col='Replicate', added_samples=None, valid_values=2):
        data = standards.merge(df, on='Accession', how='left')
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

        R2 = r2_score(Y, Y_pred)
        
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
    
class advisor:
    
    def KMeans(sorted_data, component='PC1', n_clusters=2):
    
        model = KMeans(n_clusters=n_clusters)
    
        model.fit(sorted_data[component].values.reshape([-1, 1]) )
    
        all_predictions = model.predict(sorted_data[component].values.reshape([-1, 1]))
        sorted_data['cluster'] = pd.Series(all_predictions)
        
        group = all_predictions[0]
        suggestions = sorted_data[sorted_data.cluster == group]['Sample']
        
        return sorted_data, suggestions
    
    def pca(df, standards, lfq_method='iBAQ'):
        
        stan_list = standards['Accession'].unique()  
    
        idstd = df[df.Accession.isin(stan_list)].dropna()
        
        pivot_data = idstd.dropna(subset=lfq_method).pivot_table(index='Sample', 
                                                         columns='Accession', 
                                                         values=lfq_method).fillna(0).reset_index()
        features = pivot_data.columns[1:]
        # Separating out the features
        x = pivot_data.loc[:, features].values
        # Separating out the target
        #y = pivot_data.loc[:,['Sample']].values
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
    
    def spits(df, id_col, lfq_method, replicate_dict, formatting='auto', 
             transformation=np.log2, normalization=None, valid_values=0.7, imputation='', **imp_kwargs):
        """
        """
        df.columns = df.columns.str.replace('.: ', '')
        samples = list(replicate_dict.keys())
        #conditions = list(tools.invert_dict(replicate_dict).keys())
    
        if 'Accession' not in df.columns:
            df.columns = [re.sub(id_col, 'Accession', id_col) 
                          if id_col == col else col for col in df.columns]
    
        # Removal of contaminants, and decoys
        #potential_cols = ['identified by site', 'contaminant', 'Reverse']
        #cont_key = [col for col in df.select_dtypes(exclude=np.number).columns for item in potential_cols if item in col]
    
        # Select columns with additional information 
        all_ids, priority = ['Accession'],  ['name', 'kDa']
        
        ids = [col for col in df.select_dtypes(exclude=np.number).columns if col in all_ids]
        
    
        df = df[ids+samples].replace(0, np.nan)
        log_samples = [i for i in samples if 'LOG' not in i.upper()]
        df[log_samples] = df[log_samples].apply(lambda x: transformation(x)) 
        df[samples] = clean.filter_rows_by_missing_values(df[samples], replicate_dict, valid_values)
        df = df.dropna(subset=samples, thresh=1)
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

    def mini_spits(df, std, it, id_col, added_samples=None, 
                   norm_dict= {'None':None, 'Relative': Normalization.Relative, 
                               'Median': Normalization.Median, 'Quantile': Normalization.Quantile}, 
                   transformation=np.log2, valid_values=2, std_amount_col='log2_Amount_fmol'):
        
        ids = [id_col]
        sample_cols = [i for v in it.values() for i in v]
        df = df[ids+sample_cols].replace(0, np.nan)
        log_samples = [i for i in sample_cols if 'LOG' not in i.upper()]
        df[log_samples] = df[log_samples].apply(lambda x: transformation(x)) 
        
        df = df.rename(columns={id_col:'Accession'})
        temp = df.copy()
        lfq_cols = {}
        for lfq_method in it.keys():
            try:    
                rep_dict, condition, it_2 = detective.alpacaWatson(df, it, id_col, lfq_method)
                if added_samples != None:
                    samples = [k for k, v in rep_dict.items() if v in added_samples]
                else:
                    samples = list(rep_dict.keys())
                lfq_cols[lfq_method] = samples
            except:
                continue
        filtering = ['Accession'] + [i for k, v in lfq_cols.items() for i in v]
        
        temp = temp[filtering]
        
        if norm_dict != None:
            new_cols = {}
            export = pd.DataFrame()
            for lfq_method, samples in lfq_cols.items():
                for norm_meth in norm_dict:
                    norm_temp = df[['Accession'] + samples].set_index('Accession').copy()
                    test = Normalization.normalizer(norm_temp, samples, norm_dict[norm_meth])
                    norm_lfq_method = f'{norm_meth} {lfq_method}'
                    new_samples = [i.replace(lfq_method, norm_lfq_method) for i in samples]                    
                    test.columns = new_samples
                    export = pd.concat([export, test], axis=1)
                    new_cols[norm_lfq_method] = new_samples
            lfq_cols = new_cols.copy()
            temp = export.copy()
        results = std.merge(temp, on='Accession', how='right').dropna(subset=std_amount_col)

        res_hm = pd.DataFrame([[meth, advisor.Regressor(results, cols, std_amount_col)] 
                  for meth, cols in lfq_cols.items()], columns=['Method', 'score'])
        
        res_hm[['normalization', 'method']] = res_hm['Method'].str.split(' ', expand=True, n=1)
        met = res_hm['score'].idxmax()#.iloc[0]
        lfq_advise, norm_advise = res_hm.iloc[met]['method'], res_hm.iloc[met]['normalization']
        
        return res_hm, lfq_advise, norm_advise
        
    def X_Y(df, value_col, amount_col):
        return df[amount_col].values.reshape(-1, 1), df[value_col].values.reshape(-1, 1) 
        
    def Regressor(df, value_col, amount_col, valid_values=2):
        
        data = df[['Accession', amount_col] + value_col].dropna(subset=value_col, thresh=valid_values)
        data['mean'] = data[value_col].mean(axis=1)
        X, Y = advisor.X_Y(data, 'mean', amount_col)

        linear_regressor = LinearRegression().fit(X, Y)  # create object for the class & perform linear regression
        Y_pred = linear_regressor.predict(X)  # make predictions
        R2 = r2_score(Y, Y_pred)
    
        return np.round(R2, 3)

    def caged(raw_df, standards_clean, it, id_col, lfq_menu, norm_dict, added, 
              value_filter=2, imputation='', amount=6, concentration=0.5, 
              cleaning=True, formatting=True):
        advise = []
        
        for int_meth in lfq_menu:
            rep_dict, cond, it_2 = detective.alpacaWatson(raw_df, it, id_col, int_meth)
            inv_test = {v: k for k, v in rep_dict.items()}
            rep_test = [inv_test[i] for i in added]
            for norm_meth in norm_dict.keys():
                std_test = standards_clean[['Accession', 'log2_Amount_fmol']].copy()
                
                test = advisor.spits(raw_df, id_col, int_meth, rep_dict,
                                      cleaning=cleaning, formatting=formatting, normalization=norm_dict[norm_meth],
                                      valid_values=value_filter, imputation=imputation,)
        
                mock, mock, coef, inter, R2 = Quantification.census(test, std_test, concentration=concentration, 
                                                                           in_sample=amount, lfq_col=int_meth, filter_col='Sample',
                                                                           added_samples=rep_test, valid_values=2,
                                                                           save=False)  
                        
                advise.append([int_meth, norm_meth, R2])
                
        res_hm = pd.DataFrame(advise, columns=['method', 'normalization', 'score'])#.set_index('method')
        
        met = res_hm['score'].idxmax()#.iloc[0]
        lfq_advise, norm_advise = res_hm.iloc[met]['method'], res_hm.iloc[met]['normalization']
        
        return lfq_advise, norm_advise , res_hm
    