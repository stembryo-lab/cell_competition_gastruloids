# **RNA-Seq Analysis Pipeline**

This repository contains the scripts used to preprocess, align, and perform downstream analysis of the bulk RNA-Seq dataset from the study:

*“Temporal constraints for developmental cell competition in mosaic gastruloids.”*

---

It is organized into two main folders, each representing a key stage of the RNA-Seq analysis workflow.

## Folder Structure

### `preprocessing_and_alignment/`

This folder contains all scripts and resources used for **preprocessing raw sequencing files** and **mapping reads** to a reference genome. The main goal of this step is to generate **count data** for each sample, which will later be used in downstream analyses.

Typical steps included in this folder:
- Quality control of raw FASTQ files
- Adapter trimming and filtering
- Alignment to a reference genome
- Generation of read count matrices

### `analysis/`

This folder contains **Jupyter notebooks** (or R notebooks) for **downstream analysis** of the count data generated in the preprocessing step.

Common analyses performed here:
- Data normalization and transformation
- Quality control and exploratory data analysis
- Differential gene expression analysis
- Functional enrichment analysis
- Visualization of results

## Getting Started

1. Start with the `preprocessing_and_alignment/` folder to process raw sequencing data and generate counts.
2. Move to the `analysis/` folder to explore the data, run statistical tests, and interpret biological findings.

---

Feel free to adapt this structure and extend the pipeline depending on the specific needs of your project.
