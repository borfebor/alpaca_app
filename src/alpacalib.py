#needed packages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import bokeh
from PIL import Image
from bokeh.io import output_notebook
from bokeh.io import push_notebook, show, output_notebook
from bokeh.palettes import Category20c
from bokeh.io import export_png
from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.transform import cumsum
from math import pi
output_notebook()
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns

class alpaca:
    task = ['Your wish is my command.',
            "I'll give my best on it.",
            'Brr... Brr...', 'Working that proteins, mate.',
            'Doing some crazy science!',
            "I've gathered the best team for our duty",
            'Have you ever wondered how proteins smell?',
            "It could take a bit, meanwhile why don't you listen to Bad Bunny?"]

    done = ["Got it!",
            "Done. :)",
            "Ole, ole, ole, Cholo Simeone.",
            "Easy!",
            "Looking forward for the next collab.",
            "OH, YEAH!",
            "Hi pretty, do you wanna hang out after work? ;)",
            "For Gondor!!!"]

    def importer(self):  # Opens a bar that ask you which document (tsv only) you want to import, and stores it as a DataFrame
        df = pd.read_csv(input(), sep='\t')

        df = df[df['Reverse'].isna()]
        df = df[df['Only identified by site'].isna()]
        df = df[df['Potential contaminant'].isna()]
        df = df.drop(columns=['Reverse', 'Only identified by site', 'Potential contaminant'])
        df = df.rename(columns={'Protein IDs': 'accession'})

        # iBAQ = df.iloc[:, 47:51]  # df.iloc[:, 0:11]
        # ids = df.iloc[:, 0:11]
        # iBAQ = pd.concat([iBAQ, ids], axis=1)

        return df, iBAQ

    def mq_import(
            file):  # Imports and arranges the DataFrame generated from MaxQuant/Perseus data, so it can be easily manipulated.
        df = pd.read_csv(file, sep='\t')
        df = df.rename(columns={'T: Majority protein IDs': 'Protein ID'})
        return df

    def ratios(df):  # Needed for DiaCys Calculation of REDOX profile of the sample, used in diacys()
        df1 = pd.DataFrame()
        df1['cuenta'] = df['oxidation ratio'].value_counts(bins=10, sort=False)
        df1['ocps'] = (df1['cuenta'] / df1['cuenta'].sum()) * 100
        df1['axis'] = range(5, 105, 10)
        return df1

    # Builds the whole ID Proteins dataframe merging and cleaning the triplicates. It only shows the total ID Proteins
    def idproteins(df1, df2, df3):
        i = randint(0, 7)
        print("\033[1m" + task[i] + "\033[0;0m")
        df = pd.concat([df1, df2, df3])
        df = df.drop_duplicates(subset=['Protein ID'])
        print("\033[1m" + done[i] + "\033[0;0m")
        return df

        # Builds the whole triplicates matrix

    def triplicates(df1, df2, df3):
        i = randint(0, 7)
        print("\033[1m" + task[i] + "\033[0;0m")
        df = pd.merge(df1, df2, how="left", on="Protein ID", suffixes=("_01", "_02"))
        df = pd.merge(df, df3, how="left", on="Protein ID", suffixes=("", "_03"))
        df = df.rename(columns={'Protein ID_01': 'Protein ID'})
        print("\033[1m" + done[i] + "\033[0;0m")
        return df

    # Prediction database import. Fitted for importing GP4 database from Excel.
    def predictor(file):
        assistant = ['Your wish is my command.', "I'll give my best on it.", 'Brr... Brr...',
                     'Working that proteins, mate.', 'Doing some crazy science!',
                     "I've gathered the best team for our duty",
                     'Have you ever wondered how proteins smell?',
                     "It could take a bit, meanwhile why don't you listen to Bad Bunny?"]
        i = randint(0, 7)
        df = pd.read_excel(file)
        df[['sp', 'Protein ID', 'Entry Name']] = df['Protein_ID'].str.split('|', expand=True)
        df = df.drop(columns='sp')
        df.head()
        return df

    # Match any database with results.
    def match(df1, df2):
        i = randint(0, 7)
        print("\033[1m" + task[i] + "\033[0;0m")
        df = df1[df1['Protein ID'].isin(df2['Protein ID'])]
        print("\033[1m" + done[i] + "\033[0;0m")
        return df

    # Bokeh interactive graph to check the proportions of your subcellular locations.
    def locations(df):
        gp4 = predictor('201118_GP4_SCL_complete.xlsx')

        loc = match(gp4, df)

        i = randint(0, 7)
        x = loc['Predicted_SCL'].astype(str).value_counts()
        data = pd.Series(x).reset_index(name='value').rename(columns={'index': 'location'})
        print(data)
        data['angle'] = data['value'] / data['value'].sum() * 2 * pi
        data['color'] = Category20c[len(x)]
        z = 110 * (data['value'] / data['value'].sum())
        data['value'] = z
        p = figure(plot_height=350, title="", toolbar_location=None,
                   tools="hover", tooltips="@location: @value{0.2f} %", x_range=(-.5, .5))
        p.annular_wedge(x=0, y=1, inner_radius=0.15, outer_radius=0.25, direction="anticlock",
                        start_angle=cumsum('angle', include_zero=True), end_angle=cumsum('angle'),
                        line_color="white", fill_color='color', legend='location', source=data)
        p.axis.axis_label = None
        p.axis.visible = False
        p.grid.grid_line_color = None
        show(p)

    # Takes UPS2 standards and draws the standard curve for the carculations, give an excel file with the identified UPS2 proteins. Returns 2 DataFrames(UPS2 Standards & Data with the measured UPS2 standards)

    def ups(df):
        ups2 = pd.read_excel('210205_UPS2_standards.xlsx')
        ups2['log2 Amount'] = np.log2(ups2['Amount (fmoles)'])
        ups2['Mass fraction (fmol/µg extract)'] = ups2['Amount (fmoles)'] / 10.6
        ups2['log2 Mass fract'] = np.log2(ups2['Mass fraction (fmol/µg extract)'])
        print("\033[1m" + 'How much volume - in µl - have you used to resuspend your UPS?' + "\033[0;0m")
        volume = float(input()) * 1000000  # input in µL
        ups2['Concentration (ng/µl)'] = ups2['Amount (fmoles)'] * ups2['Average MW (Da) (calculated)'] / volume
        print("\033[1m" + 'How much have you spiked in your samples?' + "\033[0;0m")
        dilution = float(input()) * ups2['Concentration (ng/µl)']
        ups2['Final concentration (ng/µl)'] = ups2['Concentration (ng/µl)'] * dilution
        # print("UPS standards information has been stored in ups2 dataframe")

        df = df.rename(columns={'Protein ID': 'accession'})
        data = pd.merge(ups2, df, on='accession', how='right')
        data = data.dropna(subset=['Amount (fmoles)'])
        X = data["log2 Amount"].values.reshape(-1, 1)  # values converts it into a numpy array
        Y = data["iBAQ"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
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
        print('Coefficient of determination: %.2f'
              % r2_score(Y, Y_pred))
        print(' ')

        sns.set(style='darkgrid', font_scale=1.5)
        g = sns.lmplot(x="log2 Amount", y="iBAQ", data=data, aspect=1.5)
        g.set_axis_labels("UPS2 amount (log2)", "UPS2 IBAQ (log2)")
        g.savefig('210310_UPS2_Quantification.png', bbox_inches='tight', pad_inches=0.5)
        with pd.ExcelWriter('ups2_quantification.xlsx') as writer:
            ups2.to_excel(writer, sheet_name='UPS standards')
            data.to_excel(writer, sheet_name='Measured UPS2')
        return ups2, data, coef, inter

    # Function used for the calculation of fmol amounts from our sample

    def quant(self):
        return lambda x: 2 ** ((x - inter) / coef)

    # Enrichment standard quantification from iBAQ quantification on MaxQuant data.

    # DETAILS
    # It takes our results DataFrame and the Excel file in which we have the Enrichment Standards.
    # In the Standards, it needs an 'Accession' column for matching results and 'fmol(Log2)' for plotting and calculating the standard curve.
    # It returns 2 DataFrame, the Spiked Standards & the Measured Standards. Both are also exported in an excel file in their respective excel sheet.

    def std_amount(df):
        e1 = df
        print("\033[1m" + 'Did you diluted your standard mix? Which dilution did you perform? (1/x)' + "\033[0;0m")
        dilution = float(input())
        print("\033[1m" + 'How much (µl) did you added into your sample? (e.g. 5.2)' + "\033[0;0m")
        volume = float(input())
        print("\033[1m" + 'How much is the initial sample volume (ml)? (e.g. 45 ml)' + "\033[0;0m")
        sample_vol = float(input())
        e1['Added Standard proteins (ng)'] = e1[['Mix concentration (µg/µl)']] / dilution * volume
        e1['Amount (fmol) in sample'] = e1['Added Standard proteins (ng)'] / e1['MW (kDa)'] * 1000
        e1['Standard fmol/ml sample'] = e1['Amount (fmol) in sample'] / sample_vol
        return e1

    def enrichment(df, file):
        std = pd.read_excel(file)

        df = df.rename(columns={'Protein ID': 'Accession'})
        data = pd.merge(df, std, on='Accession', how='right')
        data = data.dropna(subset=['fmol (Log2)'])
        X = data["fmol (Log2)"].values.reshape(-1, 1)  # values converts it into a numpy array
        Y = data["iBAQ"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
        # linear_regressor = LinearRegression().fit(X, Y)  # create object for the class & perform linear regression
        # Y_pred = linear_regressor.predict(X)  # make predictions

        # The coefficients
        # print('Coefficients: \n', linear_regressor.coef_)
        # The mean squared error
        # print('Mean squared error: %.2f'
        # % mean_squared_error(Y, Y_pred))
        # The coefficient of determination: 1 is perfect prediction
        # print('Coefficient of determination: %.2f'
        # % r2_score(Y, Y_pred))

        sns.set(style='darkgrid', font_scale=1.5)
        g = sns.lmplot(x="fmol (Log2)", y="iBAQ", data=data, aspect=1.5)
        g.set_axis_labels("STD amount (log2)", "STD IBAQ (log2)")
        g.savefig('210418_STD_Quantification.png')
        # print(linregress(x,y))
        with pd.ExcelWriter('STD_quantification.xlsx') as writer:
            std.to_excel(writer, sheet_name='Enrichment standards')
            data.to_excel(writer, sheet_name='Measured eSTD')
        return std, data