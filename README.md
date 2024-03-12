![ALPACA](https://github.com/borfebor/alpaca_app/blob/62b6915c377ccc9af4bc85aed6e715ea27c583b3/ALPACA_LOGO2.png)

> An app by [B. Ferrero-Bordera](https://www.linkedin.com/in/borjaferrero/)
### Requirements

Required packages for a proper function of the pipeline and last tested versions:

- pandas 2.1.4
- numpy 1.24.3
- matplotlib 3.8.1
- scipy 1.11.3
- scikit-learn 1.4.1.post1
- sklearn 1.0.1
- seaborn 0.11.2
- thefuzz 0.20.0


**Other packages**
- math / pi
- random / seed, randint

## Input Data Requirements

- **Protein Groups** file from MaxQuant as `.csv`, `.txt` or `.xlsx`.
- **Enrichment Standards** file as `.csv`, `.txt` or `.xlsx`. It requires 3 columns for a proper execution (Accession, MW (kDa) & Mix concentration (µg/µl). See Enrichment for more details.

#### Labwork input
- **prep**: list with our enrichment addition steps using the following order [`Dilution_1`, `Added volume_1`, `Sample Volume_1`, `Dilution_2`, `Added volume_2`, `Sample Volume_2` ... `Dilution_N`, `Added volume_N`, `Sample Volume_N`] in which N stands for the total amount of different Enrichments that we have performed with our samples.
- **enrichment_type_dict**: dictionary of lists containing which experimentals conditions have been enriched in the same way during our sample preparation. `enrichment_type_dict` = {Enrichment_1 : [Control, Treatment_1], Enrichment_2 : [Treatment_2]}
- **count_dict**: dictionary in which keys stands for our counted Experimental condition and values for the number of cells / ml that has been counted. `count_dict` = {Control : 1000000, Treatment : 2000000}

#### Optional input
- **condition**: list with the experimental conditions assayed [`Condition_1`, `Condition_2`,...,`Condition_N`] with N being the total number of experimental conditions assayed. It is created by `alpaca.experimenter()`. 
- **n**: Integer with the number of experimental conditions. It is in the output of `alpaca.experimenter()`.
- **r**: Integer with the number of replicates. Returned by `alpaca.experimenter()`.
- **rep**: List of replicate names. [`Replicate_1`, `Replicate_2`,...,`Replicate_R`] with R being the total amount of replicates. It is returned by `alpaca.replicator()` and used as an argument `alpaca.condition_format()`. 
