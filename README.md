# RNA-seq Analysis Snakemake Pipeline

An end-to-end **RNA-seq differential expression analysis pipeline** implemented using **Snakemake**, covering quality control, read alignment, quantification, statistical analysis, and functional enrichment. This workflow is designed for reproducible, modular, and scalable RNA-seq analyses.

---

## Table of Contents
- [Overview](#overview)  
- [Pipeline Features](#pipeline-features)  
- [Folder Structure](#folder-structure)  
- [Installation](#installation)  
- [Usage](#usage)  
- [Outputs](#outputs)  
- [Dependencies / Tools](#dependencies--tools)  
- [Author](#author)  

---

## Overview
This pipeline automates the entire RNA-seq workflow, starting from raw FASTQ files to downstream analysis and visualization of differential gene expression. It is suitable for **single or multiple conditions**, and is built for reproducibility using Snakemake’s workflow management.

---

## Pipeline Features
- Quality control of raw reads (e.g., FastQC)  
- Read alignment to reference genome (e.g., HISAT2/STAR)  
- Gene quantification (featureCounts / htseq-count)  
- Differential expression analysis using **DESeq2**  
- Visualization:  
  - Volcano plots  
  - MA plots  
  - PCA plots  
  - Heatmaps  
- Functional enrichment analysis (GO BP, MF, CC using clusterProfiler)  
- Export of all results in organized directories

---

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
