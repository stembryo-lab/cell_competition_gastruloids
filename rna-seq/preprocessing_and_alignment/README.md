# **RNA-Seq Analysis Workflow**

This repository contains the scripts used to preprocess, align, and perform downstream analysis of the bulk RNA-Seq dataset from the study:

*“Temporal constraints for developmental cell competition in mosaic gastruloids.”*

---

## Preprocessing & Alignment

### 1. Quality control & adapter trimming
To assess raw FASTQ quality and perform adapter trimming, run:

```bash
./qc_trimming.sh
```

### 2. Alignment with STAR
Trimmed reads are then aligned to reference genome version GRCm38.p6 using STAR aligner by running:
```bash
./mapping.sh
```

### 3. Constructing expression matrix
Once the reads are properly mapped, 'countdata_construct.ipynb' concatenates all aligned for each sample into single count matrix.

## Downstream analysis
