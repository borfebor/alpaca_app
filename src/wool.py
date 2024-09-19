#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:22:19 2024

@author: borfebor
"""

import pandas as pd

class yarnScissors:

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