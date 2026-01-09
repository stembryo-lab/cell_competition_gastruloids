# Image analysis

Image analysis for this work was primarily focused on using cellular segmentation to extract spatial information in gastruloids at the single-cell level. For this purpose, we used the software qlivecell, available at:
https://github.com/dsb-lab/qlivecell

In particular, we used version v0.8. To install this version, run:

python -m pip install 'qlivecell @ git+https://github.com/dsb-lab/qlivecell@0.8'

The qlivecell package provides tools for segmentation, error correction, and downstream analysis. Cellular segmentation relies on StarDist (https://github.com/stardist/stardist) and Cellpose (https://github.com/MouseLand/cellpose). Please refer to their respective GitHub pages for installation instructions.

Inside the cell_competition_gastruloids/image-analysis folder you will find two subfolders:

- examples
- complete_analysis

------------------------------------------------------------

examples folder

The examples folder contains a single runnable Jupyter notebook, example.ipynb, which provides a minimal, end-to-end demonstration of all image-analysis steps performed in this work using a single gastruloid.

The notebook is fully executable and is intended to illustrate the complete analysis workflow, from raw image data to quantitative, single-cell–resolved spatial descriptors.

------------------------------------------------------------

complete_analysis folder

The complete_analysis folder contains all the code used to generate the image-analysis results presented in the article.
This folder is provided for reference only, as the scripts will not run directly due to the absence of the full dataset in this repository.

Below, you will find:
1. A complete guide on how to run the example notebook.
2. A walkthrough of the structure of the complete_analysis folder, allowing the reader to map each script to its corresponding figure or panel in the paper.

------------------------------------------------------------

Example Workflow (Runnable)

The example workflow is fully contained in the notebook:
- examples/example.ipynb

This notebook provides a minimal, end-to-end workflow illustrating how 3D gastruloid images are processed, segmented, and analyzed using qlivecell. All analyses are performed on a single example gastruloid and are presented in the recommended execution order.

The example dataset consists of a single file:
- gastruloid.tif: A 3D multichannel image of a gastruloid used as a proof of concept for the full analysis pipeline.

Before running the notebook, make sure that qlivecell is installed and that the path to the example gastruloid is correctly specified.

The notebook covers the following analysis steps:

1. Cell segmentation
Initial 3D cell segmentation of the gastruloid using StarDist and Cellpose through the qlivecell interface, producing labeled segmentation masks and associated metadata.

2. Debris size threshold estimation
Estimation of a size-based threshold to distinguish real cells from small segmented objects corresponding to noise or debris.

3. Removal of segmented debris
Cleaning of the segmentation masks by removing objects below the estimated size threshold.

4. Cell-type counts
Quantification of the number of cells in each population within the gastruloid, assuming the presence of multiple cell types.

5. Nuclear fluorescence quantification
Single-cell nuclear fluorescence quantification using the segmentation masks, including background subtraction and extraction of per-cell intensity values.

6. Local cell density computation
Computation of local cell density from the 3D segmentation results, capturing spatial crowding and tissue organization.

7. Neighborhood composition analysis
Extraction of the cellular neighborhood composition for each cell, quantifying the identities of nearby cells. This analysis is only meaningful when more than one cell population is present.

8. Radial distribution analysis
Analysis of the radial distribution of cells relative to the gastruloid center, enabling the study of spatial patterning along the radial axis.

Together, these analyses illustrate how raw 3D image data are transformed into quantitative, single-cell–resolved spatial descriptors of gastruloid organization.

------------------------------------------------------------

Full Analysis Pipeline (Reference Only)

Within the image-analysis framework, each type of downstream analysis based on the segmentation results is organized into dedicated folders:

- image-analysis/code/segmentation
- image-analysis/code/apoptosis_segmentation
- image-analysis/code/counts
- image-analysis/code/density
- image-analysis/code/neighborhood_composition
- image-analysis/code/radial_distribution
- image-analysis/code/p53
- image-analysis/code/p53_degron
- image-analysis/code/YAP

For additional details, each file contains a short description at the beginning explaining its purpose.