# Docker Guide for RNA-seq Snakemake Pipeline

This guide explains how to use the provided `Dockerfile` to run the RNA-seq Snakemake pipeline in a fully reproducible container environment.

---

## Overview

The Docker image contains all necessary tools and dependencies to execute the pipeline:

- **Snakemake** – workflow manager
- **FastQC** – quality control
- **STAR** – read alignment
- **featureCounts** – gene quantification
- **R** with packages: DESeq2, clusterProfiler, EnhancedVolcano, pheatmap, ggplot2, AnnotationDbi, org.Sc.sgd.db

This allows anyone to run the pipeline without installing dependencies on their host machine.

---

## Prerequisites

- Docker installed: [https://www.docker.com/get-started](https://www.docker.com/get-started)
- Git (optional, to clone the repository)

---

## Steps to Run

## 1. Clone the repository
```bash
git clone https://github.com/Abishag24/rna-seq-analysis-snakemake-pipeline.git
cd rna-seq-analysis-snakemake-pipeline

2. Build the Docker image
docker build -t rna_seq_pipeline 

3. Run the pipeline in the container
docker run --rm -v $(pwd):/workspace -w /workspace rna_seq_pipeline snakemake --cores 4 --use-conda

## Input Requirements

Place input data in the data/ folder:

Raw FASTQ files

Reference genome (FASTA)

Gene annotation (GTF)

Modify config.yaml if you want to change paths or parameters.

## Outputs

results/ – DESeq2 results, plots, and enrichment analysis

counts/ – gene-level counts

logs/ – workflow logs

## Notes

This Docker container ensures reproducibility; the same results can be generated regardless of host system.

If the pipeline is updated, rebuild the image to include new dependencies or workflow changes:

docker build --no-cache -t rna_seq_pipeline .
