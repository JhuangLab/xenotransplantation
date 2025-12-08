# Integrated Multi-Omics Analysis of xenotransplantation 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/) [![Nextflow](https://img.shields.io/badge/nextflow-%3E%2023.10-099acd.svg)](https://www.nextflow.io) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32-brightgreen.svg)](https://snakemake.readthedocs.io)

**Repository**: https://github.com/jinyancool/MultipleMyeloma\
**Principal Investigator**: Jhuanglab\
**Contact**: hiekeen \[at\] gmail.com

## Project Overview

**MultipleMyeloma** is a comprehensive, fully reproducible research compendium supporting the study:

> **"Integrated genomic and transcriptomic multi-omics defines molecular subtypes, clonal evolution trajectories, and therapeutic vulnerabilities in xenotransplantation"**

xenotransplantation (MM) remains an incurable plasma cell malignancy despite recent therapeutic advances. Bulk sequencing approaches mask extensive intra-tumoral clonal heterogeneity and fail to capture immune microenvironment dynamics. This project addresses these limitations by performing deep **multi-omics profiling** on a large cohort of newly diagnosed and relapsed/refractory MM patients, integrating the following platforms on matched samples:

-   Whole-genome sequencing (WGS, tumor + germline)
-   Bulk RNA sequencing (polyA + total RNA)
-   10x Genomics 3' and 5' single-cell RNA sequencing (scRNA-seq)
-   10x Genomics V(D)J single-cell BCR sequencing (scBCR-seq)
-   10x Genomics V(D)J single-cell TCR sequencing (scTCR-seq)

### Key Scientific Aims

1.  Define robust molecular subtypes of MM using integrated genomic and transcriptomic features\
2.  Reconstruct clonal evolution and infer phylogenetic trees from WGS + scRNA-seq trajectories\
3.  Characterize plasma cell phenotypic states (cycling, stress response, pre-plasma, drug-resistant)\
4.  Map the bone marrow immune microenvironment and T-cell exhaustion trajectories\
5.  Reconstruct clonal and receptor-specific B-cell and T-cell dynamics using paired scBCR/TCR data\
6.  Identify actionable therapeutic vulnerabilities (e.g., synthetic lethality, immunotherapy targets)\
7.  Develop multi-omics-based prognostic and predictive models for clinical outcome and treatment response

# Quick Start

## 1. Clone the repository

git clone https://github.com/jinyancool/MultipleMyeloma.git cd MultipleMyeloma

## 2. Install conda environments

conda env create -f envs/environment_python.yml conda env create -f envs/environment_r.yml conda activate mm_multiomics

## 3. Example: Run the full Nextflow single-cell pipeline

nextflow run workflows/nf-10x/main.nf -profile conda,singularity --outdir results/scRNA

## Software Environment

-   Scanpy, Scirpy, MuData, scvi-tools, CellRank
-   Seurat, SeuratObject, Signac, ArchR
-   CellPhoneDB v4, LIANA+, NicheNet, Tensor-cell2cell
-   bcftools, GATK4, Delly, Manta, GRIDSS (WGS)
-   STAR, Salmon, kallisto, FusionCatcher (RNA-seq)
-   TensorFlow / PyTorch for deep learning-based subtype classification and trajectory inference

Detailed environments are in `envs/`.

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Funded by \[Your Funding Agency/Grant Number\]. We thank 10x Genomics, NanoString, and Akoya Biosciences for technical support and early-access programs.
