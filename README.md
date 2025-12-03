# Cellular Competition in Gastruloids

This is the repository associated with the article **"Temporal constraints for developmental cell competition in mosaic gastruloids"**, submitted to *Nature Cell Biology*. The analysis—and therefore the repository—is divided into two parts:

- **RNA sequencing analysis**  
  Correspondence: Pau Pascual Mas — EMAIL

- **Image analysis**  
  Correspondence: Pablo Casaní Galdón — pablocasani@gmail.com

## RNA sequencing analysis

TBD

## Image analysis

Image analysis for this work was mainly focused on using cellular segmentation as a way to extract spatial information in gastruloids at the single-cell level. For that, we used the software **qlivecell**, located at:  
https://github.com/dsb-lab/qlivecell

In particular, we used version `v0.8`. To install that version, run:

    python -m pip install 'qlivecell @ git+https://github.com/dsb-lab/qlivecell@0.8'

The `qlivecell` package serves as a tool to handle segmentation, error correction, and downstream analysis. Cellular segmentation relies on **StarDist** (https://github.com/stardist/stardist) and **Cellpose** (https://github.com/MouseLand/cellpose). Please refer to their respective GitHub pages for installation instructions.

Inside the `cell_competition_gastruloids/image-analysis` folder you will find two folders:

- `examples`  
- `complete_analysis`

#### `examples` folder

The `examples` folder contains minimal examples of all analyses performed in this work, using a single gastruloid in each case. These scripts are designed to be executable and illustrate the workflow step by step.

#### `complete_analysis` folder

The `complete_analysis` folder contains all the code used to produce the image-analysis results shown in the article.  
This folder is provided **for reference only**, as the scripts will not run directly because the complete dataset is not included in this repository.

Below, you will find:

1. A complete guide on how to run the examples.  
2. A walkthrough of the structure of the `complete_analysis` folder so the reader can map each script to its corresponding figure or panel in the paper.

---

### Example Workflows (Runnable)

---

### Full Analysis Pipeline (Reference Only)
