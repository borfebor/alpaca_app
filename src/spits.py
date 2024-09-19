#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 12:33:23 2024

@author: borfebor
"""

from thefuzz import fuzz
import pandas as pd
import numpy as np
import os
import re
from itertools import permutations
from scipy.stats import truncnorm
from fancyimpute import KNN, IterativeSVD, IterativeImputer, MatrixFactorization, NuclearNormMinimization, SoftImpute, SimpleFill

class detective:
    
    def alpacaHolmes(df):
        """
        Indentifies the column containing our accessions, intensity methods available and respective columns
        """
        id_col = tools.find_protein_id(df)
        is_pivot = tools.check_format(df, id_col)
        #int_columns = tools.identify_candidates(df, id_col=id_col)
        it = tools.got_it(df)
        return id_col, is_pivot, it
    
    def alpacaWatson(df, it, id_col, lfq_method):
        """
        Groups the intensity columns based on their conditions and replicates.
        """
        it = tools.got_it(df)
        it = tools.needle(it)
        candidates = tools.in_polish_please(it, lfq_method)
        export = tools.in_the_darkness_bind_them(candidates)#, test
        replicate_dict, conditions = tools.wrap_it_up(export, it, lfq_method)
        it[lfq_method] = [i for i in it[lfq_method] if i in list(replicate_dict.keys())]

        return replicate_dict, conditions, it
            #except:
            #    st.error('An error happened here')

class Normalization:

    def Quantile(df):
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
        
    def Median(series):
        return series - series.median()
    
    def Relative(series):
        pro = series.apply(lambda x: 2**x)
        return np.log2(pro / pro.sum())
        
    def normalizer(data, samples, method=None):
        if method != None:
            data[samples] = method(data[samples])
            return data
        else:
            return data

class clean:
        
    # Function to filter rows based on the missing value threshold
    def filter_rows_by_missing_values(data, replicate_dict, threshold=0.7):
        thresh = 1 - threshold 
        condition_groups = tools.invert_dict(replicate_dict)
        # Create a mask for rows to keep
        keep_rows_mask = pd.Series([True] * data.shape[0], index=data.index)
        
        for condition, columns in condition_groups.items():
            # Subset the data based on the condition columns
            subset = data[columns]
            # Calculate the fraction of missing values per row
            missing_fraction = subset.isnull().mean(axis=1)
            # Update the mask to exclude rows with missing values above the threshold
            keep_rows_mask &= (missing_fraction <= thresh)
        
        # Filter the original data based on the mask
        filtered_data = data[keep_rows_mask]
        return filtered_data

class Imputation:
    
    # Function to impute data by condition
    def by_condition(data, replicate_dict, impute_function, **kwargs):
        imputed_data = pd.DataFrame(index=data.index)
        condition_groups = tools.invert_dict(replicate_dict) 
        for condition, columns in condition_groups.items():
            # Subset the data based on the condition columns
            subset = data[columns]
            # Apply the imputation function
            imputed_subset = impute_function(subset, **kwargs)
            # Combine the imputed subset with the rest of the data
            imputed_data = pd.concat([imputed_data, imputed_subset], axis=1)
        
        return imputed_data
    
    # Function for Lowest of Detection (LOD)
    def impute_lod(data, lod=0.01):
        return data.fillna(lod)
    
    # Function for Random Drawing from a Left-Censored Normal Distribution (ND)
    def impute_nd(data, lod=0.01):
        mean, std = data.mean().mean(), data.std().mean()
        def draw_value():
            return truncnorm.rvs(a=(lod - mean) / std, b=np.inf, loc=mean, scale=std)
        
        return data.applymap(lambda x: draw_value() if pd.isnull(x) else x)
    
    # Function for k-Nearest Neighbors (kNN)
    def impute_knn(data, n_neighbors=5):
        imputer = KNN(k=n_neighbors)
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)
    
    # Function for Local Least Squares (LLS) using IterativeImputer
    def impute_lls(data, max_iter=10):
        imputer = IterativeImputer(max_iter=max_iter)
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)
    
    # Function for Random Forest (RF) using IterativeImputer and MatrixFactorization
    def impute_rf(data, max_iter=10):
        imputer = IterativeImputer(estimator=MatrixFactorization(), max_iter=max_iter)
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)
    
    # Function for Singular Value Decomposition (SVD)
    def impute_svd(data, n_components=2):
        imputer = IterativeSVD(rank=n_components)
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)
    
    # Function for Bayesian Principal Component Analysis (BPCA) using MatrixFactorization
    def impute_bpca(data, n_components=2, max_iter=100):
        imputer = MatrixFactorization(rank=n_components, max_iters=max_iter)
        return pd.DataFrame(imputer.fit_transform(data), columns=data.columns, index=data.index)

class tools:
    
    def invert_dict(replicate_dict):
        inv_map = {}
        for k, v in replicate_dict.items():
            inv_map[v.split(';')[0]] = inv_map.get(v.split(';')[0], []) + [k]
        return inv_map
    
    def check_format(df, id_col='Accession'):
        return df[id_col].value_counts().mean() == 1
    
    def find_protein_id(data):
        df = data.copy().dropna(axis=1)
        id_col_candidates = [i for i in df.select_dtypes(exclude=[np.number, 'bool']).columns if 'Protein' in i]
        id_sorter = [df[i].apply(lambda x: len(str(x.split(';')[0]))).mean() for i in id_col_candidates]
        return id_col_candidates[id_sorter.index(np.min(id_sorter))]

    def spectro_shorter(df):
        
        look_for = ['PG.QUANTITY','PG.LOG2QUANTITY', 'LOG2QUANTITY', 'LOG2', ]
        cols = []
        for p in look_for:
            c = [i for i in df.columns if p in i.upper()]
            #print(p, c)
            if len(c) > len(cols):
                cols = c
        return cols
    
    def identify_candidates(df, max_i = 1e4, id_col='Accession'):
        
        sort = df.select_dtypes(np.number).replace(0, np.nan)
        spectro_cols = tools.spectro_shorter(sort)
        
        if len(spectro_cols) > 0:
            return spectro_cols
        
        else:
            item = sort.sum() / df[id_col].nunique()
            try:
                #thresh = sort.Intensity.min()
                print(f'Thresh set to {thresh}')
            except:
                thresh = max_i
    
            challengers = list(item[item > thresh].index)
    
            candidates = [i for i in challengers if 'PG' in i]
    
            if len(candidates) > 0:
                return candidates
            elif len(challengers) > 0:
                return challengers
            
    def classify_intensity_column(column_name):
        # Define keywords and patterns for intensity methods
        intensity_keywords = ['iBAQ', 'LFQ', 'APEX', 'Quantity', 'Intensity', "Top3"]
        intensity_patterns = [r'iBAQ', "LFQ", r"APEX", r"Quantity", r"Intensity", r"Top3"]  # Regular expressions for pattern matching
        
        # Check for direct keyword matches
        for keyword in intensity_keywords:
            if keyword in column_name:
                return keyword
        
        # Check for pattern matches using regular expressions
        for pattern in intensity_patterns:
            if re.search(pattern, column_name, re.IGNORECASE):
                return re.findall(pattern, column_name, re.IGNORECASE)[0]
        
        # If no match found, return None
        return None
    
    def got_it(column_names):
        intensity_columns = {}
    
        for column_name in column_names:
            intensity_method = tools.classify_intensity_column(column_name)
            if intensity_method:
                intensity_columns.setdefault(intensity_method, []).append(column_name)
           
        return intensity_columns

    def in_polish_please(it, lfq_method):
        """
        Cleans the naming of the intensity columns
        """
        arranged, candidates = tools.racoon(it, lfq_method)
        candidates = tools.gatekeeper(candidates)
        return [f'Batman {i}' if len(i) <= 6 else i for i in candidates]
    
    def needle(it):
        exit_dict={}
        
        for lfq_method in it:
            it_clean = [re.sub(lfq_method, '', i) for i in it[lfq_method]]
            candi_sorted = sorted(it_clean, key=len)
            suff = []
            
            try:
                for i in candi_sorted:
                    diff = [re.sub(i, '', x) for x in candi_sorted if (i in x)]
                    [suff.append(i) for i in diff if i not in suff]
            except:
                return it
                
            new_intensities = list(dict.fromkeys([f'{p}{lfq_method}' if i.find(lfq_method) > i.find(p) else f'{lfq_method}{p}' for i in it[lfq_method] for p in suff ]))
            new_intensities = [i for i in new_intensities 
                        if tools.has_numbers(re.sub('Log2', '', i,  flags=re.IGNORECASE)) == False]
            if len(new_intensities) != len(candi_sorted):
                to_arrange = it[lfq_method].copy()
            
                for i in sorted(new_intensities, key=len, reverse=True):
                    fount = [col for col in to_arrange if i in col]
                    exit_dict[i] = fount
                    [to_arrange.remove(value) for value in fount]
    
        return exit_dict
    
    def in_the_darkness_bind_them(candidates):
        """
        Groups the intensity columns by condition by their likelyhood to belong to the same condition
        """
        test = pd.DataFrame([[x, y, fuzz.ratio(x, y), [x,y], len(x), fuzz.ratio(re.sub(r'\d', '', x), re.sub(r'\d', '', y))] 
                             for x, y in permutations(candidates,2) if fuzz.ratio(x, y) > 85], columns=['x', 'y', 'ratio', 'combined', 'len_x', 'ratio_no_digits'])
                
        test['common'] = test['combined'].apply(lambda x: len(os.path.commonprefix(x)))
        test['shared'] = test.common / test.len_x
        export = []
        for i in test.x.unique():
            #sorted_i = test[(test.x == i) & (test.len_x - test.common < (test.len_x*0.15))]
            sorted_i = test[(test.x == i) & (test.shared >= 0.85)]
            mx = sorted_i.ratio_no_digits.max()# - sorted_i.ratio.std()
            the_good_ones = list(sorted_i[sorted_i.ratio_no_digits == mx].y.unique())
            the_good_ones.append(i)
            the_good_ones.sort()
            the_good_ones = [re.sub('Batman', '', i).strip() for i in the_good_ones]
            if the_good_ones not in export:
                export.append(the_good_ones)
        return export
    
    def wrap_it_up(export, it, lfq_method):
        """
        Groups the intensity columns by condition.
        Returns a dictionary and a list.
            - Dictionary summarizes sample columns, conditions and replicates (e.g., {'LFQ intensity Control_1': ['Control', 'Replicate_1']}, 
            - The list contains the identified conditions (e.g., ['Control', 'Treatment'])
        """  
        sample_dict = {}
        conditions = []
        for i in export:
            #common = os.path.commonprefix(i)
            common = re.sub(r'_0$', '', os.path.commonprefix(i))
            conditions.append(common.strip('_ '))
            samples = [[col, f"{common.strip('_ ')};{sample.replace(common, '').strip('_ ')}"] 
                   for sample in i for col in it[lfq_method] if sample in col]
            sample_dict.update(samples)
        return sample_dict, conditions
    
    def remove_dates(candidates):
        return [re.sub(r'^\d{6}_', '', i)
                 for i in candidates]
        
    def clear_the_path(columns):
        return [col.split('\\')[-1] for col in columns]
    
    def remove_numbering(columns):
        return [re.sub(r'\[\d+\]', '', i).strip(' ') for i in columns]
    
    def remove_substrings(s):
        s = s.replace(".", "")
        s = s.replace("PG", "")
        s = s.replace("raw", "")
        return s
    
    def remove_the_tail(columns):
        return [tools.remove_substrings(i).strip(' ') for i in columns]
    
    def has_numbers(inputString):
        return bool(re.search(r'\d', inputString))
    
    def racoon(it, lfq_method):
        lfq = [i for i in it[lfq_method]# if lfq_method in i 
               if tools.has_numbers(re.sub('Log2', '', i,  flags=re.IGNORECASE)) == True]
        
        common_prefix = os.path.commonprefix(lfq)
        
        if lfq_method in common_prefix:
            lfq = [re.sub(common_prefix, '', s).strip() 
                   if len(re.sub(common_prefix, '', s).strip()) > 1
                   else re.sub(' '.join(common_prefix.split(' ')[:-1]), '', s).strip()
                   for s in lfq]
            #print(common_prefix)
            return common_prefix, lfq
        
        reversed_strings = [s[::-1] for s in lfq]
        common_suffix = os.path.commonprefix(reversed_strings)[::-1]
        
        if lfq_method in common_suffix:
            lfq = [re.sub(common_suffix, '', s).strip() 
                   if len(re.sub(common_suffix, '', s).strip()) > 1
                   else re.sub(' '.join(common_suffix.split(' ')[1:]), '', s).strip()
                   for s in lfq]
            #print(common_suffix)
            return common_suffix, lfq
    
    def gatekeeper(columns):
        return tools.remove_the_tail(tools.remove_dates(tools.clear_the_path(tools.remove_numbering(columns))))