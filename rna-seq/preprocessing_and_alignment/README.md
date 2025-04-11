# Preprocessing and Alignment

This folder contains scripts and a notebook for **preprocessing raw RNA-Seq data**, performing **read alignment** to a pre-built genome index, and **constructing a count matrix** for downstream analysis.

## Contents

### `qc_trimming.sh`
This script performs:
- **Quality control** with **FastQC**
- **Adapter and quality trimming** using **Trim Galore**

Run this script on raw FASTQ files to generate cleaned reads for alignment.

---

### `mapping.sh`
This script runs **STAR** to:
- **Align trimmed reads** to a **pre-built reference genome index**
- **Generate gene-level counts** using the `--quantMode GeneCounts` option

> ⚠️ Assumes that the STAR genome index is already built and the path is correctly specified in the script.

---

### `countdata_construct.ipynb`
This notebook:
- Loads individual `ReadsPerGene.out.tab` files from STAR
- Extracts count data for each sample
- Concatenates the results into a single **gene-by-sample count matrix**
- Outputs a clean matrix ready for downstream analysis

Includes comments and modular sections for clarity and reproducibility.

## Usage Workflow

1. **Run `qc_trimming.sh`** to prepare trimmed reads.
2. **Run `mapping.sh`** on the trimmed reads to obtain per-sample gene counts.
3. **Execute `countdata_construct.ipynb`** to merge individual count files into a single matrix.

---

### Requirements
- Tools: `FastQC`, `Trim Galore`, `STAR`
- Input: Raw FASTQ files, STAR index (pre-built), and genome annotation (e.g. GTF file)

Make sure paths in the scripts are set correctly for your working environment.
