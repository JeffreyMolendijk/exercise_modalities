## Introduction

Code repository for the manuscript: 

<cite>Phosphoproteomics of three exercise modalities identifies canonical exercise signaling and C18ORF25 as an AMPK substrate regulating skeletal muscle function.</cite>

<br>

## Usage

* Clone this repository to your local computer (`git clone https://github.com/JeffreyMolendijk/exercise_modalities.git`)
* Install the R-packages as used in the individual scripts. 
* Execute the scripts `human_phosphosite.R`, `human_protein.R`, `c18orf25_phosphosite.R` and `c18orf25_protein.R`.
* Inspect the output tables and images in `data/export/`

Each of the named folders contains the scripts required to replicate different sections of the analysis (human phosphoproteomics, human proteomics, c18orf25 knockout phosphoproteomics and c18orf25 knockout proteomics). To locate the scripts of interest, please refer to the folder structure below.

<br>

## Folder structure

```
exercise_modalities
│   README.md    
│
└───data/
│   │
│   └───input/
│   │   │   c18orf25_phospho.txt
│   │   │   c18orf25_protein.txt
│   │   │   human_SPS_v2.RData
│   │   │   Phospho (STY)Sites_Filter.xlsx
│   │   │   Physiological individual data - teto-USE THIS.xlsx
│   │   │   proteinGroups_filterRename.txt
│   │
│   └───export/
│       │
│       └───human_phosphosite/
│       │   ...
│       │   
│       └───human_protein/
│       │   ...
│       │
│       └───c18orf25_phosphosite/
│       │   ...
│       │
│       └───c18orf25_protein/
│       │   ...
│
└───R/
    │
    └───human_phosphosite/
    │   │   human_phosphosite.R
    │   │   ...
    │
    └───human_protein/
    │   │   human_protein.R
    │   │   ...
    │
    └───c18orf25_phosphosite/
    │   │   c18orf25_phosphosite.R
    │   │   ...
    │
    └───c18orf25_protein/
        │   c18orf25_protein.R
        │   ...
```

<br>

## Input data

[INSERT DESCRIPTION OF ALL INPUT TABLES > ANY PRE-PROCESSING?]

| filename                  | description                               |
| -------------             | -------------                             |
| c18orf25_phospho.txt      | c18orf25 knockdown phosphoproteomics data |
| c18orf25_protein.txt      | c18orf25 knockdown proteomics data        |
| human_SPS_v2.RData        | object containing human phosphosites      |

<br>

## Analysis outputs

### `human_phosphosite`
* PCA plot (before / after processing)
* Kinase enrichment analysis plots and data
* Limma DE results (phospho_grand_DE.txt)
* Limma adjusted p-value distribution
* FunScor annotated table (site_annotations.csv)
* Phosphosite-trait correlations (phosphosite_trait_correlation_spearman.csv)
* Processed data table (human_phosphosite_rba.csv)

### `human_protein`
* PCA plot (before / after processing)
* Limma DE results (prot_grand_DE.txt)
* Limma adjusted p-value distribution
* Processed data table (human_protein_rba.csv)

### `c18orf25_phosphosite`
* Limma DE results (phospho_grand_DE.txt)

### `c18orf25_protein`
* Limma DE results (prot_grand_DE.txt)
* GSEA result tables


<br>

## Citation

<cite>TO BE ADDED</cite>

<br>

## Contact
For more information, please contact Benjamin L. Parker (myemail@unimelb.edu.au).