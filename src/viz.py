#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 21:00:40 2023

@author: borfebor
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from random import randint
from math import pi
import seaborn as sns
from src.alpaca_new import alpaca
from st_aggrid import AgGrid
import altair as alt
import time

class Viz:
    
    def boxplot(df, categorical, numerical, color):
        
        chart = alt.Chart(df).mark_boxplot(extent='min-max').encode(
            x=categorical,
            y=numerical,
            color=color
            )
        
        return chart
    
    def identifications(df, categorical, color, grouper):
        
        df_grouped = df.groupby(grouper)['Accession'].nunique().reset_index()
        
        
        bar = alt.Chart(df_grouped).mark_bar().encode(
            y = alt.Y(categorical),
            x = alt.X('Accession:Q'),
            color = color
            )
        
        text = alt.Chart(df_grouped).mark_text(dx=-15, dy=3, color='white').encode(
                x=alt.X('Accession:Q'),
                y=alt.Y(categorical),
                #detail='site:N',
                text=alt.Text('Accession:Q', format='.0f')
            )
        
        chart = bar + text
        
        return chart
    
    
    def Regression(df, amount, lfq_method, R):
        
        fitting = f'R2 = {round(R, 3)}'
            
        reg = alt.Chart(df).mark_point().encode(
                alt.X(amount,
                      scale=alt.Scale(domain=(df[amount].min()-df[amount].std()/2,
                                              df[amount].max()+df[amount].std()/2))),
                alt.Y(lfq_method,
                      scale=alt.Scale(domain=(df[lfq_method].min()-df[lfq_method].std()/2, 
                                              df[lfq_method].max()+df[lfq_method].std()/2))),
                ).properties(title=fitting, height=500)
        
        
        chart = reg + reg.transform_regression(amount, lfq_method).mark_line() 
        
        return chart
    
    def displot(df, lfq_method):
    
        chart = alt.Chart(df).mark_bar(
            opacity=0.3,
            binSpacing=0
        ).encode(
            alt.X(lfq_method, bin=alt.Bin(maxbins=100)),
            alt.Y('count()', stack=None),
            alt.Color('Sample:N')
            )
        
        return chart
    
    def scatter(df, x, y, variance, columns):
        
        name = [f'PC{var[0]+1} ({100*round(var[1],2)} %)' for var in enumerate(variance)]
        
        pos_x = columns.index(x)
        pos_y = columns.index(y)       
    
        chart = alt.Chart(df).mark_circle(size=60).encode(
            x=alt.X(x, axis=alt.Axis(title=name[pos_x])),
            y=alt.Y(y, axis=alt.Axis(title=name[pos_y])),
            color='Condition',
            tooltip='Sample'
            ).interactive()
        
        return chart
    
    def z_score(data, intensity_method='LFQ'):

        data['z_score'] = np.nan
    
        for sample in data.Sample.unique():
    
            value = data[data.Sample == sample][intensity_method]
            mean = value.mean()
            sd = value.std()
            data['z_score'] = np.where(data.Sample == sample, 
                                       (data[intensity_method] - mean)/sd, 
                                        data.z_score)
            
        sns.catplot(data, x='Condition', y='z_score', kind='box', hue='Replicate')
        return data
    
    def heatmap(df, x, y, c, z_score=False, color_scheme='redblue'):
        
        source = df.copy().dropna(subset=y)
        
        if z_score == True:
        
            source = Viz.z_score(source, c)
            c = 'z_score'
        
        quant_75 = source[c].quantile(0.75)
        quant_50 = source[c].quantile(0.5)
        quant_25 = source[c].quantile(0.25)
        
        chart = alt.Chart(source).mark_rect().encode(
            alt.X(x),
            alt.Y(y),
            alt.Color(c, scale=alt.Scale(
            domain=[quant_25,quant_50,quant_75], 
            scheme=color_scheme, 
            #interpolate=method
            ),
        legend=alt.Legend(direction='horizontal', orient='top', title=None)
        )
        )
        return chart
