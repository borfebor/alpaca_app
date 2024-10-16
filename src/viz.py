
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 21:00:40 2023

@author: borfebor
"""

import numpy as np
import seaborn as sns
import plotly.graph_objs as go
import plotly.express as px
from sklearn.linear_model import LinearRegression

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
                    xaxis_title=None,
                    hovermode="x",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.2,
                        xanchor="right",
                        x=1
                ))  
            chart.update_yaxes(automargin=True)
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
                    xaxis_title=None,
                    hovermode="x",
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.2,
                        xanchor="right",
                        x=1
                )) 
            box.update_yaxes(automargin=True)
        except:
            box = f'Oops! Apparently I cannot plot based on {categorical}'
        return box
    
    
    def Regression(df, amount, lfq_method, R):
        
        # regression
        reg = LinearRegression().fit(np.vstack(df[amount]), df[lfq_method].values)
        df['bestfit'] = reg.predict(np.vstack(df[amount]))
        
        fitting = f'R² = {round(R, 3)}'
        try: 
            # plotly figure setup
            chart=go.Figure()
            chart.add_trace(go.Scatter(name='Standard proteins', x=df[amount], y=df[lfq_method].values, mode='markers'))
            chart.add_trace(go.Scatter(name='line of best fit', x=df[amount], y=df['bestfit'], mode='lines', line_color='Black'))
            
            # plotly figure layout
            chart.update_layout(xaxis_title = 'fmol of standard (log2)', yaxis_title = f"{lfq_method} (log2)")
            chart.update_traces(marker=dict(
                                color='LightSkyBlue',
                                size=8,
                                line=dict(width=2,
                                        color='DarkSlateGrey')),
                      selector=dict(mode='markers'))
            chart.update_layout(title_text=fitting, title_x=0.5,  width=500,
             height=600)
            chart.update_yaxes(automargin=True) 

        except:
            chart = "Oops! Something went wrong when I tried to plot your standards :("
                  
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
            chart.update_yaxes(automargin=True)
        except:
            chart = 'Sorry :( This is a bit embarrassing but something went wrong'
        
        return chart
    
    def scatter(df, x, y, variance, columns):
        try: 
            name = [f'PC{var[0]+1} ({100*round(var[1],2)} %)' for var in enumerate(variance)]
    
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
                    xaxis_title=name[int(x[-1])-1],
                    yaxis_title=name[int(y[-1])-1],
                    legend=dict(
                        orientation="h",
                        yanchor="bottom",
                        y=1.02,
                        xanchor="right",
                        x=1
                )) 
            chart.update_yaxes(automargin=True)
        except:
            chart = 'Sorry :( This is a bit embarrassing but something went wrong'
        
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
            chart.update_layout(
                    xaxis_title=None,
                    yaxis_title=None) 
            chart.update_yaxes(automargin=True)
        except:
            chart = 'Sorry :( This is a bit embarrassing but something went wrong'

        return chart


