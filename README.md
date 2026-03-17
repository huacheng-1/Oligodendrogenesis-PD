# Oligodendrogenesis-PD
## Contents

- [Overview](#overview)
- [Repository Structure](#Repository-Structure)
- [System Requirements](#system-requirements)
- [Installation Guide](#Installation-Guide)
- [Data](#Data)
- [License](#License)


# Overview

This repository contains the R analysis scripts associated with the manuscript: Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease.

# Repository Structure

├─[LDSC](./figure1a&b_ldsc): `figure1 a&b` code.

├─[EWCE](./figure1c_EWCE):: `figure1 c` code.

├─[scRNA-seq](./figure2scRNA-seq_analysis):`figure2` code.

├─[Bulk RNA-seq](./figure3a-c&e-g_Bulk_RNA-seq):`figure3 a-c&e-g` code.

├─[AUCell](./figure3d&h_AUCell):`figure3 d&h` code.

└─[OL subclusters](figure3i-l_OL-subclusters):`figure3 i-l` code.

```
├── figure1a&b_ldsc/
│   ├── fast_match2_2.pl                              # Match SNPs against cell-type annotation files for LDSC input
│   ├── step1_bed_generation/
│   │   ├── integrated_scRNA_to_bed.R                 # Convert scRNA-seq cell-type markers to BED files for LDSC
│   │   └── scATAC_to_bed.R                           # Extract scATAC-seq peaks per cell type & export as BED files for LDSC
│   ├── step2_ldsc.txt                                # Bash pipeline: merge baseline annotations, generate .annot files, compute LD scores & run stratified LDSC
│   └── step3_plot.R                                  # Compute p-values from LDSC z-scores & visualize cell-type enrichment
├── figure1c_EWCE/
│   └── figure1c(EWCE).R                              # EWCE cell-type enrichment test for PD DEGs
├── figure2scRNA-seq_analysis/
│   ├── figure2a-c/
│   │   ├── GSE243639_pre-process.R                   # QC, normalization, clustering & annotation for GSE243639 (human PD)
│   │   ├── GSE178265_pre-process.R                   # QC, normalization, clustering & annotation for GSE178265 (human PD)
│   │   ├── GSE193688_pre-process.R                   # QC, normalization, clustering & annotation for GSE193688 (human PD)
│   │   ├── figure2a.R                                # Integrates GSE243639/178265/193688 with Harmony and unified cell-type annotation
│   │   └── figure2b,c.R                              # UMAP, proportion bar chart, marker dot plot & expression heatmap from integrated object
│   ├── figure2d-f.R                                  # snRNA-seq full pipeline (QC → integration → clustering → annotation → visualization) for CNP0000892 α-synuclein PD mouse snRNA-seq
│   └── figure2g-i.R                                  # snRNA-seq full pipeline (QC → integration → clustering → annotation → visualization) for MPTP mouse snRNA-seq
│
├── figure3a-c&e-g_Bulk_RNA-seq/
│   ├── microarray_data_pre-process/
│   │   ├── GSE7621.R                                 # Download GSE7621 & run limma DEG (PD vs CON)
│   │   └── GSE49036.R                                # Download GSE49036 & run limma DEG (PD vs CON)
│   ├── figure3a_e_volcano.R                          # Volcano plot of limma DEGs
│   ├── figure3b_f_GO.R                               # GO-BP enrichment on DEGs
│   └── figure3c_g_heatmap.R                          # Heatmap of target gene expression
│
├── figure3d&h_AUCell/
│   └── figure3d&h_AUCell.R                           #Score oxidative stress gene set activity per cell with AUCell 
│
└── figure3i-l OL subclusters/
    ├── figure3i-j_OL subclusters.R                   # Sub-cluster oligodendrocytes, compare PD/CON subcluster proportions
    └── figure3k-l_GO.R                               # GO-BP enrichment on subcluster 3 and 6 
```

# System Requirements
## Hardware requirements

Oligodendrogenesis-PD requires only a bioinformatics workstation with sufficient RAM to support high-throughput in-memory genomic computations.

CPU: 13th Gen Intel(R) Core(TM) i7-13700 (2.10 GHz)

RAM: 64.0 GB


## OS Requirements
This package is supported for Windows and Linux. The package has been tested on the following systems:

Windows: Windows 11(25H2)

Linux: Ubuntu 22.04

## Software Dependencies
### Required R Packages/Linux software
Bioconductor version 3.19,

DESeq2 version 1.44.0,

enrichplot version 1.24.4,

clusterProfiler version 4.12.6,

msigdbr version 24.1.0,

fgsea version 1.30.0,

limma version 3.60.6,

Signac version 1.14.9002,

SeuratObject version 5.3.0,

harmony version 1.2.3,

samtools version 1.21,

cellranger version 8.0.0,

FastQC version 0.12.1,

hisat2 version 2.2.1,

LDSC version 1.0.1 (https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format ).

# Installation Guide
R version 4.4.1 installation guide:
1.InstallR from https://cran.r-project.org/
2.(0ptional)Install RStudio from https://posit.co/download/rstudio-desktop/
3.Install all required packages by running:
```
# Install dependencies
# Install R from https://cran.r-project.org/
# Typical installation time: Approximately 10-20 minutes on a standard desktopcomputer with a stable internet connect(longer if Seurat dependencies # areLnot pre-installed).
#**From CRAN:**
install.packages(c(
"tidyverse", # version >= 1.3.0
"ggplot2"，# version >= 3.4.0
"dplyr", # version >= 1.1.0
"readr", # version >= 2.1.0
"patchwork", # version >= 1.1.0(for multi-panel figures)
"RColorBrewer", # version >= 1.1.3
"Bioconductor", # version >= 3.19
"DESeq2", # version >= 1.44.0
"enrichplot", # version >= 1.24.4
"clusterProfiler", # version >= 4.12.6
"msigdbr", # version >= 24.1.0
"fgsea", # version >= 1.30.0
"limma", # version >= 3.60.6
"Signac", # version >= 1.14.9002
"SeuratObject", # versioin >= 5.3.0
"harmony" # version >= 1.2.3
))
#**From Bioconductor (single-cell analysis):**
if (Irequire("BiocManager", quietly = TRUE))install.packages("BiocManager")

```


# Data
Four PD GWAS summary statistics datasets derived primarily from individuals of European ancestry were compiled. These included: Leal et al., 2025 (doi: 10.1101/2025.07.18.25331793) ; GP2 Parkinson’s Disease GWAS, European ancestry, 2025; Kim et al., 2024 (doi: 10.1038/s41588-023-01584-8);  and Nalls et al., 2019 (doi: 10.1016/S1474-4422(19)30320-5).

6 human PD datasets (GSE243639, GSE178265, GSE193688, GSE147672, GSE7621 and GSE49036) were obtained from the GEO database. The MPTP mouse dataset was retrieved from Figshare (identifier: https://figshare.com/articles/dataset/Single cell_sequencing_screening_of_biomarkers_associated_with_Lipophagy_in_the_MPTP_model_mice_of_Parkinson_s_disease/26724526/1?file=48569005), and the α-synuclein A53T transgenic mouse dataset was obtained from the CNGBdb (CNSA; accession code CNP0000892, https://db.cngb.org/cnsa/).

# License
This code is released under the **MIT License**. See LICENSE' file for details.

