# Downstream Analysis

This folder contains notebooks used for the **downstream analysis** of RNA-Seq data. These analyses are performed using the count data generated in the `preprocessing_and_alignment/` step.

## Overview

Each notebook in this folder is focused on a specific analysis task and can be **executed independently**. Notebooks are well-documented with **descriptions and comments** to guide the user through each step of the analysis process.

### Typical Analyses Included:
- **Data exploration and quality control**
- **Normalization and filtering**
- **Principal component analysis (PCA) and clustering**
- **Differential gene expression (DGE) analysis**
- **Integration to embryo psedobulk**
- **Spatial reconstruction of cell populations**
- **Visualization of results**

## Structure

Each notebook typically follows this structure:
1. **Load libraries and data**
2. **Preprocess input (e.g., filtering low counts)**
3. **Apply specific analysis (e.g., DGE, PCA)**
4. **Generate and interpret plots**
5. **Save results if applicable**

## How to Use

- You can run each notebook independently depending on the analysis you wish to perform.
- Make sure required input files (e.g., count matrix, metadata) are available or adjust paths as necessary.
- Follow comments within each notebook for a better understanding of methods and assumptions used.

---

If you have any questions or need help with a specific notebook, feel free to open an issue or reach out.

