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

import plotly.graph_objs as go
import plotly.io as pio
from plotly.subplots import make_subplots
import seaborn as sns
import plotly.express as px
import statsmodels.api as sm


class Viz:
    
    def identifications(df, categorical, color, grouper):
        
        colores = df[color].unique()
            
        colors = dict(zip(colores, sns.color_palette('Set1', len(colores)).as_hex()))
        try: 
            df_grouped = df.groupby(grouper)['Accession'].nunique().reset_index()
  
            chart = px.bar(df_grouped, 
                           y='Accession', x=categorical,
                          color=color,
                          text_auto=True,
                          #line_color='#000000',
                          color_discrete_sequence=px.colors.qualitative.Set3,
                          ).update_traces(marker_line_width=2)
            chart.update_layout(
                    yaxis_title='ID proteins',
                    hovermode="x",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.2,
                        xanchor="right",
                        x=1
                ))  
        except:
            error_on = (',').join(group for group in grouper)
            chart = f'Oops! Apparently I cannot plot based on the given categories ({error_on})'
                
        return chart

    def boxplot(df, categorical, numerical, color, lfq_method='Intensity'):
        try:
            colores = df[color].unique()
    
            categories = df[categorical].unique()
            
            colors = dict(zip(colores, sns.color_palette('Set3', len(colores)).as_hex()))
            
            box = go.Figure()
            previous = ''
                    
            for num, group in enumerate(categories):
                
                color = [(colors[key], key) for key in colors if key in group][0]
    
                if color[1] == previous:
                    legend=False
                else:
                    legend=True
                
                print(group)
                previous = color[1]
                box.add_trace(go.Box(
                            y=df[df[categorical] == group][numerical],
                            x=[group] * len(df[df[categorical] == group][numerical]),
                            name=color[1],
                            fillcolor=color[0],
                            line_color='#000000',
                            showlegend=legend
                        ))
      
            box.update_layout(
                    yaxis_title=lfq_method,
                    hovermode="x",
                )
                
            box.update_layout(
                    yaxis_title=numerical,
                    hovermode="x",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.2,
                        xanchor="right",
                        x=1
                )) 
        except:
            box = f'Oops! Apparently I cannot plot based on {categorical}'
        return box
    
    
    def Regression(df, amount, lfq_method, R):
        try:
            fitting = f'R² = {round(R, 3)}'
    
            chart = px.scatter(
                        x=df[amount], 
                        y=df[lfq_method],
                        trendline='ols',
                        trendline_color_override='Black',
                        hover_name=df['Accession'],
                        labels={ 
                                "x": "fmol of standard (log2)",  "y": f"{lfq_method} (log2)",
                        },
                width=500,
                height=600
            )
            chart.update_traces(marker=dict(
                                color='LightSkyBlue',
                                size=8,
                                line=dict(width=2,
                                        color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
            chart.update_layout(title_text=fitting, title_x=0.5)
        except:
            box = f'Oops! Something went wrong when I tried to plot your standards :('
        
        return chart
    
    def displot(df, lfq_method):
        try:
            chart = px.histogram(df, 
                       x=lfq_method, 
                       color="Condition",
                       hover_data=df.columns,
                       color_discrete_sequence=px.colors.qualitative.Set3, 
                       facet_col="Condition",
                       facet_row='Replicate',
                       nbins=25,
                       barmode='overlay').update_traces(marker_line_width=2)
            chart.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    
            chart.update_layout(
                    bargap=0,
                    hovermode="x",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.2,
                        xanchor="right",
                        x=1
                ))
        except:
            chart = f'Sorry :( This is a bit embarrassing but something went wrong'
        
        return chart
    
    def scatter(df, x, y, variance, columns):
        try: 
            name = [f'PC{var[0]+1} ({100*round(var[1],2)} %)' for var in enumerate(variance) if x in var[1] if y in var[1]]
    
            chart = px.scatter(df,
                        x=x, 
                        y=y,
                        color='Condition',
                        symbol="Condition",
                        color_discrete_sequence=px.colors.qualitative.Set3,
                width=500,
                height=600
                ).update_xaxes(
                    zeroline=True,
                    color='black'
                ).update_yaxes(
                    zeroline=True,
                    color='black'
                )
            
            chart.update_traces(marker=dict(
                                size=10,
                                line=dict(width=2,
                                        color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
    
            chart.update_layout(
                    xaxis_title=name[0],
                    yaxis_title=name[1],
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                )) 
        except:
            chart = f'Sorry :( This is a bit embarrassing but something went wrong'
        
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
        try:
            source = df.copy().dropna(subset=y)
            
            if z_score == True:
            
                source = Viz.z_score(source, c)
                c = 'z_score'
            
            quant_75 = source[c].quantile(0.75)
            quant_50 = source[c].quantile(0.5)
            quant_25 = source[c].quantile(0.25)
    
            hm = source.pivot_table(columns=y, index=x, values=c)
    
            chart = px.imshow(hm,
                              text_auto=True, aspect="auto",
                             color_continuous_scale=color_scheme)
        except:
            chart = f'Sorry :( This is a bit embarrassing but something went wrong'

        return chart
