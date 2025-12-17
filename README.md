# RNA-seq Analysis Snakemake Pipeline

An end-to-end **RNA-seq differential expression analysis pipeline** implemented using **Snakemake**.  
The workflow is designed to be **reproducible, modular, and scalable**, and demonstrates standard RNA-seq data processing practices.

---

## Overview
This pipeline automates RNA-seq analysis starting from raw FASTQ files through read alignment, gene quantification, and differential expression analysis.

Some downstream analyses (visualization and functional enrichment) are currently under active development.

---

## Pipeline Features
- Quality control of raw reads (FastQC)
- Read alignment to reference genome (STAR)
- Gene-level quantification (featureCounts)
- Differential expression analysis using **DESeq2** *(in progress)*
- Visualization (PCA, volcano plots, heatmaps) *(planned)*
- Functional enrichment analysis using **clusterProfiler** *(planned)*

---

## Folder Structure

```text
├── workflow/        # Snakemake workflow files (Snakefile, config)
├── scripts/         # R scripts for DESeq2 and visualization
├── data/            # Reference genome and annotation files
├── results/         # Generated outputs (ignored in Git)
├── logs/            # Log files (ignored in Git)
├── envs/            # Conda environment YAML files
└── README.md
