# oligodendrogenesis-Parkinson-s-disease
Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease
## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#Installation-Guide)
- [Instructions for Use](#Instructions-for-Use)
- [Citation Database](#Citation-Database)


# Overview

This repository provides the computational framework and analytical pipelines used in our study: “Aberrant oligodendrogenesis in the substantia nigra promotes oxidative stress-dependent PANX1 activation and neurodegeneration in Parkinson’s disease.” By integrating multi-omics analyses with cell type-specific genetic manipulation to identify OLs as critical drivers of PD pathogenesis. Partitioned heritability analysis revealed significant enrichment of PD genetic risk in OLs, and single-nucleus RNA sequencing demonstrated consistent OL expansion in the SN of human postmortem tissue and two independent mouse models. Transcriptomic integration identified OLs as the primary cellular source of ROS in the PD SN.

# Repo Contents

- [LDSC](./figure1_a&b_ldsc): `LDSC` code.
- [EWCE](./figure1_c_EWCE): `EWCE` code.
- [scRNA-seq](./figure2&figure3_i-l_scRNA-seq_analysis): `scRNA-seq` code.
- [Bulk RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq): `Bulk RNA-seq` code.
- [AUCell](./figure3_d&h_AUCell): `AUCell` code.


# System Requirements
## Hardware requirements

oligodendrogenesis-Parkinson-s-disease requires only a standard computer with enough RAM to support the in-memory operations.
CPU: 13th Gen Intel(R) Core(TM) i7-13700 (2.10 GHz)
RAM: 64.0 GB

## OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:
Windows: Windows 11(25H2)
Linux: Ubuntu 22.04

## Environment requirements
Bulk RNA-seq analysis: R version 4.4.1, Bioconductor version 3.19, DESeq2 version 1.44.0, enrichplot version 1.24.4, clusterProfiler version 4.12.6, msigdbr version 24.1.0, fgsea version 1.30.0, limma version 3.60.6.
ATAC-seq analysis: R version 4.5.1, Signac version 1.14.9002.
ScRNA-seq analysis: R version 4.5.1, SeuratObject version 5.3.0, harmony version 1.2.3, samtools version 1.21, cellranger version 8.0.0, FastQC version 0.12.1, hisat2 version 2.2.1.
Partitioned Heritability Analysis: LDSC version 1.0.1 (https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format ).

# Installation Guide
R version 4.4.1 installation guide:
```
# Install dependencies
# Ubuntu/Debian
sudo apt update
sudo apt install -y build-essential gfortran libreadline-dev xorg-dev libbz2-dev \
    liblzma-dev libpcre2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev libicu-dev

# CentOS/RHEL
sudo yum groupinstall "Development Tools"
sudo yum install -y gcc-gfortran readline-devel libX11-devel bzip2-devel \
    xz-devel pcre2-devel libcurl-devel openssl-devel zlib-devel libicu-devel

# Download R 4.4.1 source code
wget https://cran.r-project.org/src/base/R-4/R-4.4.1.tar.gz
tar -xzvf R-4.4.1.tar.gz
cd R-4.4.1

# Configuration and compilation
./configure --prefix=/usr/local/R-4.4.1 --enable-R-shlib
make -j$(nproc)  # 使用所有CPU核心加速编译
sudo make install
#--prefix=/usr/local/R-4.4.1 指定安装路径。
#--enable-R-shlib 允许 R 被其他程序调用（如 RStudio）。

# Set environment variables
echo 'export PATH=/usr/local/R-4.4.1/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# Verify installation
R --version
```

# Instructions for Use
*LDSC can run normally in this document file. If there are any errors, please check the configuration environment.
- [LDSC](./figure1_a&b_ldsc): 
  This document contains the process code for LDSC operation, including a series of codes for extracting differentially expressed genes from single-cell data, generating baseline data, calculating heritability, and p-values. There are three steps in total, step 1 and step 3 are the code that runs in the `R` environment. Step 2 is to run the LDSC package program in the Linux environment.
- [EWCE](./figure1_c_EWCE):
  This document contains scoring codes for single-cell gene sets, including data preprocessing and specificity calculation, EWCE analysis function construction, differential gene extraction, and scoring based on background gene sets.
- [scRNA-seq](./figure2&figure3_i-l_scRNA-seq_analysis):
  This document contains processing code for five sets of single-cell data, and the pre-processing of the integrated data is separately placed in[整合数据前处理](./figure2&figure3_i-l_scRNA-seq_analysis/数据整合前处理).The processing procedures for single-cell data of CNP0000892, MPTP data, and integrated data are named as follows:[CNP0000892](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_d&e&f(CNP).R), [MPTP](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_g&h&i(MPTP).R), [整合数据](./figure2&figure3_i-l_scRNA-seq_analysis/figure2_a&b&c(整合数据).R). The GO pathway enrichment analysis process for integrating data subcluster codes and subclusters is saved in directory figure3_i&j(亚聚类).R与figure3_k&l(GO).R.
- [Bulk_RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq):
  This document contains the processing code for two sets of simple RNA sequencing data, and the differential gene analysis process file is placed in[前处理](./figure3_a-c&e-g_Bulk_RNA-seq/前处理). The running process code for the GO chart of volcano map, heat map, and enrichment pathway analysis are all placed in[Bulk_RNA-seq](./figure3_a-c&e-g_Bulk_RNA-seq).
- [AUCell](./figure3_d&h_AUCell):
  This document contains scoring codes for the integrated single-cell data of two sets of simple RNA sequencing data.

# Citation Database

Four PD GWAS summary statistics datasets derived primarily from individuals of European ancestry were compiled. These included: Leal et al., 2025 (doi: 10.1101/2025.07.18.25331793) ; GP2 Parkinson’s Disease GWAS, European ancestry, 2025; Kim et al., 2024 (doi: 10.1038/s41588-023-01584-8);  and Nalls et al., 2019 (doi: 10.1016/S1474-4422(19)30320-5).
6 human PD datasets (GSE243639, GSE178265, GSE193688, GSE147672, GSE7621 and GSE49036) were obtained from the GEO database. The MPTP mouse dataset was retrieved from Figshare (identifier: https://figshare.com/articles/dataset/Single cell_sequencing_screening_of_biomarkers_associated_with_Lipophagy_in_the_MPTP_model_mice_of_Parkinson_s_disease/26724526/1?file=48569005), and the α-synuclein A53T transgenic mouse dataset was obtained from the CNGBdb (CNSA; accession code CNP0000892, https://db.cngb.org/cnsa/).
