# Kraft et. al 2022, Workflow

**Patient-specific identification of genome-wide DNA-methylation differences between intracranial and extracranial melanoma metastases**

This is the pipeline that we used for the analysis of patient-matched DNA-methylation data of melanoma metastases.

# Instructions

This repository has five parts: 

* **MainScript.R**: The main script to use is MainScript.R . It includes all analysis and can be run in R. 
* **HelperFunctions.R**: All functions that are used in the main script are stored in the separate file, HelperFunctions.R. Please load all these helper functions to your R workspace before applying the pipeline.
* **Folder: Annotations**: Contains all functional annotations that are used by the pipeline to analyze the methylation data.
* **Folder: Figures**: This folder contains all scripts to create the main figures that are part of the manuscript. A separate folder is used for each figure. Each folder contains an .RData file that comprises all data required to create the specific figure. These R.Data files are created in the pipeline (see MainScript.R). Each folder additionally contains a separate R-script that creates the figure in PDF format.
* **Folder: Programs**: This folder contains the Java JAR file implemetation of the HMM by Seifert et. al, 2014 that was used to make the predictions for the methylation states of individual CpGs (see MainScript.R).

All required data to run the pipeline is available from Zenodo under [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7074484.svg)](https://doi.org/10.5281/zenodo.7074484)

 Please download these data from Zenodo. The following data sets are used in the pipeline:
* Table_methylationSamples.csv 
* Table_availableExpressionData.tsv 


# License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
