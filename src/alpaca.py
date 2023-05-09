#Needed libraries
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bokeh

from random import seed
from random import randint

from math import pi
import altair as alt
import pandas as pd
from PIL import Image


from bokeh.io import push_notebook, show, output_notebook
from bokeh.palettes import Category20c
from bokeh.layouts import row
from bokeh.plotting import figure
from bokeh.transform import cumsum
output_notebook()

from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns


class alpaca:

    def importer(upload):
        m = upload.name
        if m.endswith('.txt'):
            df = pd.read_csv(upload, sep='\t')
        elif m.endswith('.csv'):
            df = pd.read_csv(upload, sep=',')
        elif m.endswith('.xlsx'):
            df = pd.read_excel(upload)

        return df

    def formater(df):
        columns = list()

        for x in df.columns:
            if 'N:' in x:
                new_string = x.replace('N: ', '')
                columns.append(new_string)
            elif 'T:' in x:
                new_string = x.replace('T: ', '')
                columns.append(new_string)
            else:
                columns.append(x)

        df.columns = columns

        for x in df.columns:
            if 'Protein ID' in x:
                df = df.rename(columns={x : 'Accession'})
                columns.append('Accession')
        placeholder = st.empty()
        default = list()
        # Checking for data cleaning
        if 'Reverse' in df.columns:
            default.append('Reverse')
        elif 'Only identified by site' in df.columns:
            default.append('Only identified by site')
        elif 'Potential contaminant' in df.columns:
            default.append('Potential contaminant')
            placeholder.warning("Your data could contain Possible Contaminants, Reverses and Decoys")

        if 'Accession' in df.columns:
            default.append('Accession')
        if 'Gene names' in df.columns:
            default.append('Gene names')
        if 'Unique peptides' in df.columns:
            default.append('Unique peptides')
        if 'Mol. weight [kDa]' in df.columns:
            default.append('Mol. weight [kDa]')


        for x in columns:
            if 'iBAQ' in x:
                default.append(x)

        return df, columns, default

    @st.cache
    def data_cleaner(df):

        if 'Reverse' in df.columns:
            df = df.loc[lambda df: df['Reverse'].isna()]
            df = df.drop(['Reverse'], axis=1)
        elif 'Only identified by site' in df.columns:
            df = df.loc[lambda df: df['Potential contaminant'].isna()]
            df = df.drop(['Potential contaminant'], axis=1)
        elif 'Potential contaminant' in df.columns:
            df = df.loc[lambda df: df['Only identified by site'].isna()]
            df = df.drop(['Only identified by site'], axis=1)

        return df

    @st.cache
    def log_transform(df):
        ibaq = [ x for x in df.columns if 'iBAQ' in x]
        df[ibaq] = np.log2(df[ibaq])
        return df

    def ups(df, ups2, volume, amount):

        ups2['log2 Amount'] = np.log2(ups2['Amount (fmoles)'])
        ups2['Mass fraction (fmol/µg extract)'] = ups2['Amount (fmoles)'] / 10.6
        ups2['log2 Mass fract'] = np.log2(ups2['Mass fraction (fmol/µg extract)'])
        ups2['Concentration (ng/µl)'] = ups2['Amount (fmoles)'] * ups2['Average MW (Da) (calculated)'] / (
                    volume * 1000000)

        dilution = amount * ups2['Concentration (ng/µl)']
        ups2['Final concentration (ng/µl)'] = ups2['Concentration (ng/µl)'] * dilution

        data = pd.merge(ups2, df, on='Accession', how='right')
        data = data.dropna(subset=['Amount (fmoles)'])
        st.markdown('**UPS2 Standards present in your data**')
        st.write(data)
        X = data["log2 Amount"].values.reshape(-1, 1)  # values converts it into a numpy array
        Y = data["iBAQ"].values.reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column
        linear_regressor = LinearRegression().fit(X, Y)  # create object for the class & perform linear regression
        Y_pred = linear_regressor.predict(X)  # make predictions

        col1, col2 = st.columns((0.5,1))
        # The coefficients
        coef = linear_regressor.coef_
        inter = linear_regressor.intercept_
        col1.metric('Coefficients:', round(float(coef), 3))
        col1.metric('Intercept:', round(float(inter), 3))

        # The mean squared error
        col1.metric('Mean squared error:',
                    round(float(mean_squared_error(Y, Y_pred)), 3))
        # The coefficient of determination: 1 is perfect prediction
        col1.metric('Coefficient of determination:',
                    round(float(r2_score(Y, Y_pred)), 3))

        fig = alt.Chart(data).mark_circle(size=60).encode(
                alt.X('log2 Amount', scale=alt.Scale(zero=False)),
                alt.Y('iBAQ', scale=alt.Scale(zero=False, padding=1)),
                ).properties(
            title='UPS2 calibration curve',
            width=500, height=300
                )

        final_plot = fig + fig.transform_regression('log2 Amount', 'iBAQ').mark_line()

        col2.altair_chart(final_plot, use_container_width=True)

        return ups2, data, coef, inter

    def condition_format(df, condition, rep):

        condition_col = list()
        condition_col = [ x for x in df.columns if condition in x]
        try:
            cd = df[['Accession', 'Gene names', 'Mol. weight [kDa]']]
        except:
            st.exception("Sorry, I couldn't found the columns Accession', 'Gene names', 'Mol. weight [kDa]")

        cd[rep] = df[condition_col]
        cd['Mean'] = cd[rep].mean(axis=1)
        cd['Condition'] = condition
        cd = cd.dropna(subset=rep, thresh=2)
        return cd

    def experimenter(df):
        condition = list()
        try:
            condition = [ x[5:-3] for x in df.columns if '_' in x]
        except:
            st.exception("I'm sorry, I could't predict your experimental conditions.")

        n = len(set(condition))  # defines the number of condition in the experiment
        r = int(len(condition) / n)
        condition = list(set(condition))

        return n, r, condition

    def replicator(n):
        """
        Gets the amount of replicates from experimenter() and creates tags for the columns.
        It helps solving the problem of variable amount of replicates in different experiments.
        It returns a list of replicate names with numbers, that will be used later by condition_format() function.
        """
        i = 1
        rep = list()
        for x in range(n):
            replicate = 'Replicate_' + str(i)
            rep.append(replicate)
            i += 1
        return rep

    def std_amount(df, list, i):  # i -> Enrichment types in our experiment, that's the reason that the function is in comprised in a loop
        e1 = df.copy()
        """
        This function will be comprised in the enricher() function, as a tool for calculating the enrichment amounts in
        every condition.
        """
        e1['Enrichment_type'] = 'Enrichment_' + str(i)
        e1['Added Standard proteins (ng)'] = e1[['Mix concentration (µg/µl)']] / list[0] * list[1]
        e1['Amount (fmol) in sample'] = e1['Added Standard proteins (ng)'] / e1['MW (kDa)'] * 1000
        e1['Standard fmol/ml sample'] = e1['Amount (fmol) in sample'] / list[2]
        return e1

    def enricher(std, clean_data, prep, enrichments, n, condition):
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
        std_list = list()
        y = 0
        z = 3
        x = 1
        for i in range(enrichments):
            tope = prep[y: z]
            cond = alpaca.std_amount(std, tope, x)
            std_list.append(cond)
            if x < enrichments + 1:
                x += 1
            y += 3
            z += 3
        e1 = pd.concat(std_list)
        std_condition = list()
        for p in range(n):
            cond = clean_data[clean_data['Condition'] == condition[p]]
            std_condition.append(cond)
        std = pd.concat(std_condition)
        data = pd.merge(std, e1, on='Accession', how='inner')
        return data

    def enrichment_type_dict(enrichments): #out-of-date function, this dict has been integrated in another part
        enrichment_type = list()
        i = 1
        enrichment_type = ['Enrichment_' + str(i + 1) for i in range(enrichments)]
        enrichment_type_dict = dict.fromkeys(enrichment_type)
        return enrichment_type_dict

    def enrichment_filter(std_data, enrichment_type_dict):
        filtered = list()
        for key, value in enrichment_type_dict.items():
            for v in value:
                std_data_filtered = std_data[(std_data['Enrichment_type'] == key) & (std_data['Condition'] == v)]
                filtered.append(std_data_filtered)
        std_data = pd.concat(filtered)
        return std_data

    def enrichment_factors(data_std, enrichment_type_dict):
        enrichment_factors = dict()
        for key in enrichment_type_dict.keys():
            EF = data_std.Enrichment[data_std['Enrichment_type'] == key].mean()
            enrichment_factors[key] = EF
        number = 1
        for key, value in enrichment_factors.items():
            st.markdown('Enrichment Factor for ' + key + ': **' + str(round(value, 3)) + '**')
            number += 1
        return enrichment_factors

    def sample_volumes(enrichment_type_dict, prep):
        volume_matcher = list()
        i = 2
        for key, value in enrichment_type_dict.items():
            volume_matcher.append([key, prep[i]])
            i += 3
        return volume_matcher

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
        dummy['Absolute_amount_fmol'] = dummy['Enrichment_factor'] * dummy['Amount (fmol)']
        dummy['Molecules'] = dummy['Absolute_amount_fmol'] * (
                6.023 * (10 ** 8))  # Avogadro's number fixed for fmol (-15)
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