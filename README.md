# RNA-seq Analysis Snakemake Pipeline

An end-to-end **RNA-seq differential expression analysis pipeline** implemented using **Snakemake**, covering quality control, read alignment, quantification, statistical analysis, and functional enrichment. This workflow is designed for reproducible, modular, and scalable RNA-seq analyses. Pipeline tested locally with public Saccharomyces cerevisiae RNA-seq data.

**Status:**  Tested locally with public Saccharomyces cerevisiae RNA-seq data  
**Note:** Input data must be provided by the user.

## Table of Contents
- [Overview](#overview)  
- [Pipeline Features](#pipeline-features)  
- [Folder Structure](#folder-structure)  
- [Installation](#installation)  
- [Usage](#usage)
- [Sample Plots Preview](#Sample-Plots-Preview)
- [Outputs](#outputs)  
- [Dependencies / Tools](#dependencies--tools)  
- [Author](#author)  


## Overview
This pipeline automates the entire RNA-seq workflow, starting from raw FASTQ files to downstream analysis and visualization of differential gene expression. 
Raw FASTQ
   │
   ▼
Quality Control (FastQC)
   │
   ▼
Read Alignment (STAR)
   │
   ▼
Gene Quantification (featureCounts)
   │    
   ▼
Differential Expression Analysis (DESeq2)
   │
   ▼
Visualization
 ├─ Volcano Plot
 ├─ MA Plot
 ├─ PCA Plot
 └─ Heatmaps
   │
   ▼
Functional Enrichment (GO BP, MF, CC)
   │
   ▼
Results Export (CSV/PNG)
It is suitable for **single or multiple conditions**, and is built for reproducibility using Snakemake’s workflow management. 

## Pipeline Features
- Quality control of raw reads (FastQC)  
- Read alignment to reference genome (STAR)  
- Gene quantification (featureCounts)  
- Differential expression analysis using **DESeq2**  
- Visualization:  
- Volcano plots  
- MA plots  
- PCA plots  
- Heatmaps  
- Functional enrichment analysis (GO BP, MF, CC using clusterProfiler)  
- Export of all results in organized directories


## Folder Structure
```text
├── workflow/                 # Snakemake workflow files (Snakefile, config)
├── counts/                   # Gene count files
├── results/                  # DE analysis results and plots
├── scripts/                  # R scripts (DESeq2, visualization)
├── data/                     # Reference genomes, annotation files
├── logs/                     # Snakemake log files
├── envs/                     # Conda environment YAML files
└── README.md                 # Project README
```

## Installation
- Required Tools
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Conda](https://docs.conda.io/en/latest/) 
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [STAR](https://github.com/alexdobin/STAR)
- [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
- [R](https://www.r-project.org/) with the following packages:
- DESeq2
- clusterProfiler

# Usage

1. Place your input files in the `data/` folder:
- Raw FASTQ files
- Reference genome (FASTA)
- Gene annotation (GTF)

### Option 1 – Simple (default Snakefile paths)

```bash
snakemake --cores 4 --use-conda
```

2. Option 2 : Update the `config.yaml` file with your file paths.

3. Run Snakemake:
```bash
snakemake --use-conda --configfile config.yaml --cores 4
```
### Input Data

Users must provide:
- Raw FASTQ files
- Reference genome (FASTA)
- Gene annotation (GTF)

Place these files in the `data/` folder before running the pipeline.

## Sample Plots Preview 
   Volcano Plot
   [#VolcanoPlot]
   MA Plot
   [#MAPlot] 
   PCA Plot
   [#PCAPlot]
   Heatmap (Top DEGs)
   [#Heatmap]

## Outputs
The workflow generates:

- `results/fastqc/` – FastQC HTML and zip reports
- `counts/` – gene-level count tables from featureCounts
- `results/deseq2/` – differential expression tables, plots, and summaries
- `logs/` – Snakemake execution logs

  # DESeq2 R script workflow:

  Reads raw count matrix and sample metadata
  Filters low-count genes and replaces NAs
  Performs DE analysis and calculates log2 fold changes
  Produces significant upregulated, downregulated, and combined DEGs
  Generates volcano, MA, PCA plots, and heatmaps
  Performs GO enrichment analysis (BP, MF, CC) using clusterProfiler
  Saves all results in results/ directory
## Dependencies / Tools
- Snakemake
- Conda 
- FastQC
- STAR
- featureCounts
- R (DESeq2, clusterProfiler)
- Python 3.x (for Snakemake)

## Author
Abishag Jacquline  
[GitHub](https://github.com/Abishag24) | [LinkedIn](https://www.linkedin.com/in/abishagjacquline/)

⚠️ Note: Runtime files, editor artifacts, and analysis outputs are excluded from version control to ensure full reproducibility.
