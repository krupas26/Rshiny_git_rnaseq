# shiny
This shiny dashboard is for **RNA Seq data analysis**. 
The following steps can be accomplished using this app:
1. Exploratory Data analysis / Quality Control
2. Differential Expression Analysis
3. Visualization of DEGs (MA plots, volcano plots, gene lists)
4. Functional Analysis - Gene Ontology (Biological Process, Cellular Component, & Molecular Function) and/or KEGG pathway
5. Gene Set enrichment Analysis

### Files
1. ui.R and server.R - Files for UI and server functions, respectively.
2. functions.R - This file contains some general functions used in ui.R and server.R.
3. environment.yml - This file can be used to build a conda environment containing all the packages used by this application. 