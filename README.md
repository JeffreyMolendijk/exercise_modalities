## Introduction

Code repository for the manuscript:

> <cite>[Phosphoproteomics of three exercise modalities identifies canonical signaling and C18ORF25 as an AMPK substrate regulating skeletal muscle function](https://www.sciencedirect.com/science/article/pii/S1550413122003023)</cite>

<br>

## Usage

- Clone this repository to your local computer (`git clone https://github.com/JeffreyMolendijk/exercise_modalities.git`)
- Install the R-packages as used in the individual scripts.
- Execute the scripts `human_phosphosite.R`, `human_protein.R`, `c18orf25_phosphosite.R` and `c18orf25_protein.R`.
- Inspect the output tables and images in `data/export/`

Each of the named folders contains the scripts required to replicate different sections of the analysis (human phosphoproteomics, human proteomics, c18orf25 knockout phosphoproteomics and c18orf25 knockout proteomics). To locate the scripts of interest, please refer to the folder structure below.

> Note: In the human experiments, the exercise modality _Strength_ is referred to as _Resistance_ in the manuscript.

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
│   │   │   human_phospho.xlsx
│   │   │   human_trait.xlsx
│   │   │   human_protein.txt
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

| filename             | description                                                                                |
| -------------------- | ------------------------------------------------------------------------------------------ |
| human_phospho.xlsx   | human exercise phosphoproteomic data                                                       |
| human_protein.txt    | human exercise proteomic data                                                              |
| human_trait.xlsx     | human plasma metabolites and muscle glycogen content                                       |
| c18orf25_phospho.txt | c18orf25 wild-type vs knockdown knockout with or without stimulation phosphoproteomic data |
| c18orf25_protein.txt | c18orf25 knockdown wild-type vs knockout proteomics data                                   |
| human_SPS_v2.RData   | object containing human phosphosites                                                       |

<br>

## Analysis outputs

### `human_phosphosite`

- PCA plot (before / after processing)
- Kinase enrichment analysis plots and data
- Limma DE results (phospho_grand_DE.txt)
- Limma adjusted p-value distribution
- FunScor annotated table (site_annotations.csv)
- Phosphosite-trait correlations (phosphosite_trait_correlation_spearman.csv)
- Processed data table (human_phosphosite_rba.csv)

### `human_protein`

- PCA plot (before / after processing)
- Limma DE results (prot_grand_DE.txt)
- Limma adjusted p-value distribution
- Processed data table (human_protein_rba.csv)

### `c18orf25_phosphosite`

- Limma DE results (phospho_grand_DE.txt)
- Kinase enrichment analysis plots

### `c18orf25_protein`

- Limma DE results (prot_grand_DE.txt)
- GSEA result tables

<br>

## Citation

> Ronnie Blazev, Christian S. Carl, Yaan-Kit Ng, Jeffrey Molendijk, Christian T. Voldstedlund, Yuanyuan Zhao, Di Xiao, Andrew J. Kueh, Paula M. Miotto, Vanessa R. Haynes, Justin P. Hardee, Jin D. Chung, James W. McNamara, Hongwei Qian, Paul Gregorevic, Jonathan S. Oakhill, Marco J. Herold, Thomas E. Jensen, Leszek Lisowski, Gordon S. Lynch, Garron T. Dodd, Matthew J. Watt, Pengyi Yang, Bente Kiens, Erik A. Richter, Benjamin L. Parker (2022). Phosphoproteomics of three exercise modalities identifies canonical signaling and C18ORF25 as an AMPK substrate regulating skeletal muscle function. Cell Metabolism, Volume 34, Issue 10.
> [https://doi.org/10.1016/j.cmet.2022.07.003](https://doi.org/10.1016/j.cmet.2022.07.003)

<br>

## Contact

For more information, please contact Benjamin L. Parker (ben.parker@unimelb.edu.au).
