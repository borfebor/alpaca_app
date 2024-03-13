![Instructions](https://github.com/borfebor/alpaca_app/blob/main/instructions.png)

> An app by [B. Ferrero-Bordera](https://www.linkedin.com/in/borjaferrero/)

This is a GUI from the Python package [Alpaca_proteomics](https://github.com/borfebor/alpaca_proteomics). The app aims to facilitate the analysis of proteomics samples for Absolute Proteome Quantification based on protein standards. 

# Tutorial

## Alpaca GUI tutorial
 
Alpaca is a GUI optimised for the analysis and quantification of proteomics samples using protein standards (e.g., [UPS2](https://www.sigmaaldrich.com/DE/en/product/sigma/ups2)). It uses the Python pipeline [Alpaca_proteomics](https://github.com/borfebor/alpaca_proteomics) for quantifying Mass Spectrometry (MS)-based proteomics samples.

### ![Photo 1. Home Page]

When accessing the GUI, a short description of the preferable workflow for the pipeline is shown in the main body. In the side bar, an uploader widget allows the users to import their ProteinGroups. CSV, TXT and Excel files are compatible with the software. Below the widget, different example datasets can be selected to explore the GUI capabilities.

Sample names and replicates for each studied condition are infered by the software based on the given names during the search on MaxQuant. 

### ![Photo 2. Structure of the GUI]

After uploading the ProteinGroups file, the GUI changes to the data analysis page. This page consist of two main parts, the body and the side bar.

- The **GUI side bar** controls quantification parameters as desired MS intensity method, normalization and protein standards. The side bar contains 2 main parts (`Data preprocessing` and `Quantification standards`). The desired data processing parameters can be choose in the `Data preprocessing` section. This include removal of contaminants, intensity normalization and data formatting. Additionally, `Quantification standards`section allows to select the used protein standards for quantification and the added amounts. By default it is set for UPS2, but custom standards can be added).

- The **GUI body** which is divided in 3 parts (`Experimental set-up`, `Data visualization` and `Your data`) allows to user to add details on the experiment and visualize the data. Details on the sample preparation can be added on the `Experimental set-up`tab. Summary statistics of the data and its quantification can be explored in the `Data visualization`tab. Finally, `Your data` corresponds to the quantified proteins in tabular format to ready to be exported.


# GUI side bar

### Data preprocessing

### Quantification standards

# GUI body

### Experimental set-up

This tab controls the sample preparation. It is an optional input in case the user wants to trace back to the original amount in the samples, rather than just quantifying what was measured in the Mass Spectrometer.

![Photo 3. Experimental set-up](https://github.com/borfebor/alpaca_app/blob/main/Screenshots/set_up.png)

#### Sample preparation

**Table 1.** Example of the input table to define the parameters used for sample parameters. An example table can be generated in the GUI based on the added conditions. This example table can be exported as a template and edited to be imported later into the GUI.
| Condition   | SampleVolume | ProteinConcentration | AmountMS | CellsPerML | TotalCultureVolume | ProteinSRM | fmolSRM | Enrichment | EnrichmentDirection | StdDilution | StdVolume |
|-------------|--------------|----------------------|----------|------------|--------------------|------------|---------|------------|---------------------|-------------|-----------|
| Cond1_t0    | 2.31         | 2.99                 | 9.67     | 4.54       | 7.54               | TNAMLN     | 4.44    | False      | Down                | 3.96        | 1.22      |
| Cond2_t1    | 2.50         | 0.20                 | 4.10     | 5.13       | 2.62               | AJFVYC     | 4.85    | True       | Down                | 2.43        | 1.51      |
| Cond3_t2    | 7.38         | 6.56                 | 2.77     | 3.66       | 3.80               | BYEKSC     | 9.71    | True       | Down                | 5.71        | 8.53      |

> - **Condition**: Condition in which the parameters were applied.
> - **SampleVolume**: Protein extract volume (µL) used for protein digestion.
> - **ProteinConcentration**: Determined protein concentration (µg/µl) in the sample.
> - **AmountMS**: Protein amount (µg) injected in the MS.
> - **CellsPerML**: Measured cells per mL of culture.
> - **TotalCultureVolume**: Total cultivation volume (µL).
> - **ProteinSRM** `Optional`: If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). This
> corresponds to the accession of measured protein in SRM to calculate
> the enrichment.
> - **fmolSRM** `Optional`: If the enrichment of a subcellular fraction has been calculated using targeted proteomics (SRM). Fmol of the
> proteins measured in the targeted proteomics measurements. 
> - **Enrichment** `Optional`: Boolean (True or False). Samples that have been enriched should be specified as True
> - **EnrichmentDirection** `Optional`: UP or DOWN. 
> - **StdDilution** `Optional`: This parameter specifies how many times the stock solution of enrichment standards has been diluted before
> adding it to the sample. If the standards were not diluted before
> addition, specify 1. Only used when the enrichment is calculated
> through the function alpaca.gathers() details of the preparation of
> the used proteins should be added. 
> - **StdVolume** `Optional`: Volume of enrichment standards (µL) added to the sample. Only used in case the enrichment is calculated through
> the function alpaca.gathers() details of the preparation of the used
> proteins should be added.

#### Proteome fraction enrichment (Optional)

In case the sample preparation includes the enrichment of a subproteomic fraction and protein standards have been used to quantify this enrichment step, the standards can be uploaded in this section. As in the sample preparation tab, a template can be generated and downloaded to edit. The GUI accepts the format described in table 2. 

**Table 2.** Example of the enrichment standards input.
| Accession | MW (kDa) | StdConcentration (µg/µl) |
|-----------|---------:|-------------------------:|
| P02768    |     10.1 |                     2.5  |
| Q9Y6K9    |     65.8 |                     0.8  |
| P05067    |     32.5 |                     1.2  |
| O75475    |     48.2 |                     3.0  |
| Q00653    |     20.9 |                     2.0  |

### Data visualization

![Photo 4. Visualization tab](https://github.com/borfebor/alpaca_app/blob/main/Screenshots/viz_option.png)
This tab allows to explore some summary statistics of the data. As it is the amount of proteins quantified (`Quantified proteins`) per sample and summary statistics of the measured intensities per sample represented as boxplot (`intensities`). Additionally, principal component analysis (`PCA`) and intensity distribution (`Distribution plot`) were added the assesment of the samples quality. The fitting of calibration standards can be visualized in `Calibration curve`. Finally, a `heatmap` can be drawn and exported if desired.

### Your Data

This section corresponds to the quantified data allowing the user to export it for further analysis. The data is exported as CSV. The data can be pivoted if the export format is not the desired, as by default it specified a row per protein quantified in each sample. Quantified data can also be filtered and sorted to export just a selected set of proteins. 

# Example datasets

# Cite us


