This code repository is associated to the publication of Cheng et al., RNA (2017). 
Source code of the main processing steps and to reproduce the figures of the paper is provided here.

## Usage
### Paper figure code
Code to generate paper figures are in the folder `src/Paper`. Files are named after the corresponding figure numbers in the paper. [You can also follow this link.](https://github.com/gagneurlab/Manuscript_Cheng_RNA_2017/tree/master/src/Paper)

### To run the code
1. [Download the repository from this link](https://github.com/gagneurlab/Manuscript_Cheng_RNA_2017/archive/master.zip) (~46M)
2. Unzip with `unzip <downloaded file>`
3. You need to start R in the root directory of this repository. 

## R Session Information
A sessionInfo file is provided `src/sessionInfo.txt`, which is the R environment on which this analysis was performed/tested.

library dependencies:

```
library(data.table)
library(LSD)
library(gridExtra)
library(ggplot2)
library(Biostrings)
library(grid)
library(gridBase)
library(dplyr)
library(multtest)
```
Also install the lastest `ggpval` package to automatically add P-values to ggplot.

```
# install.packages("devtools")
devtools::install_github("s6juncheng/ggpval")
```

## Key analysis steps
Key analysis step that are provided in `src/*.R` files. 

### Regression (with codons, tAI and the final joint model)
`src/Scer_regression.R`

`src/Pombe_regression.R`

`src/Scer_regression_Neymotin.R`

`src/compare_stAI_with_coef.R`

### Knockout strain analysis
`src/sun_mutation_analysis.R`

### Single-nucleotide effect analysis
1. `src/prepare_data-sc.R`
2. `src/snv-perturbation.R`
3. `src/snv-perturbation-gc-length.R`
4. `src/motif_single_base_join.R`
5. `src/analyze-snv-perturbation.R`

### qPCR motif validation analysis
`src/Paper/Figure4_motif_validation.R`

