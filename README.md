# Bulk-rnaseq-analysis

## Overview
Differential gene expression analysis pipeline using R

# Quick Start - Environment Setup

## Prerequisites
- `micromamba` installed ([docs](https://mamba.readthedocs.io))

## Setup 

```bash
# Clone repo 
git clone git@github.com:DavideCapozzi/bulk-rnaseq-analysis

# Move to env directory
cd bulk-rnaseq-analysis/env/

# Create environment
micromamba env create -n rnaseq_analysis -f environment.yml

# Activate
micromamba activate rnaseq_analysis

# Install R packages
Rscript install_r_packages.r 2>&1 | tee r_install.log
```

## Generate renv.lock (optional, for reproducibility)

```bash
R -e "renv::init(bioconductor='3.22')"
```

## Verify Installation

```bash
R -e "library(DESeq2); library(edgeR); library(fgsea); cat('✓ Ready\n')"
```

## Run Analysis

```bash
# All tools ready in one environment:
# - Bioinformatics: salmon, fastqc, samtools, bedtools
# - R packages: DESeq2, edgeR, ggplot2, fgsea, etc.
# - Python: numpy, pandas, scipy, jupyter

Rscript analysis/01_de_analysis.r
```

---

**Files included:**
- `environment.yml` → Complete conda environment
- `install_r_packages.r` → R/Bioconductor packages
- `renv.lock` → (optional) R package versions

## Methods
- Alignment-free quantification: Salmon
- DE tool: DESeq2
- Enrichment: fgsea (GSEA), clusterProfiler (ORA)
