# Image analysis

Image analysis for this work was primarily focused on using cellular segmentation to extract spatial information in gastruloids at the single-cell level.  
For this purpose, we used the software **qlivecell**, available at:  
https://github.com/dsb-lab/qlivecell

In particular, we used version `v0.8`. To install this version, run:

```bash
python -m pip install 'qlivecell @ git+https://github.com/dsb-lab/qlivecell@0.8'
```

The `qlivecell` package provides tools for segmentation, error correction, and downstream analysis.  
Cellular segmentation relies on **StarDist** (https://github.com/stardist/stardist) and **Cellpose** (https://github.com/MouseLand/cellpose). Please refer to their respective GitHub pages for installation instructions in case they are not properly installed during `qlivecell` installation.

Inside the `cell_competition_gastruloids/image-analysis` folder you will find two subfolders:

- `examples`  
- `complete_analysis`

---

#### `examples` folder

The examples folder contains three runnable Jupyter notebook, example1.ipynb, example2.ipynb, and example3.ipynb which provides a minimal, end-to-end demonstration of all image-analysis steps performed in this work using a set of three gastruloids

The notebooks are fully executable and are intended to illustrate the complete analysis workflow, from raw image data to quantitative, single-cell–resolved spatial descriptors. The contents of each of the examples is described below.

---

#### `complete_analysis` folder

The `complete_analysis` folder contains all the code used to generate the image-analysis results presented in the article.  
This folder is provided **for reference only**, as the scripts will not run directly due to the absence of the full dataset in this repository.

---

Below, you will find:

1. A complete guide on how to run the example scripts.  
2. A walkthrough of the structure of the `complete_analysis` folder, allowing the reader to map each script to its corresponding figure or panel in the paper.

---

### Example Workflows (Runnable)

We have designed a set of three runnable examples that represented all the image analysis work shown in the article associated with this repository. The contents of the three example are the following:

#### Example 1

This first example, found at `examples/example1.ipynb` in the form of a jupyter notebook, provides a minimal, end-to-end workflow illustrating how 3D gastruloid images are processed, segmented, and analyzed using qlivecell. All analyses are performed on a single example gastruloid and are presented in the recommended execution order.

The example dataset consists of a single file:
- `gastruloid1.tif`: A 3D multichannel image of a gastruloid containning two cell populations with two separate nuclear markers, mCherry and emiRFP.

Before running the notebook, make sure that  `qlivecell` is installed and that the path to the example gastruloid is correctly specified.

The notebook covers the following analysis steps:

1. **Cell segmentation.**
Initial 3D cell segmentation of the gastruloid using StarDist and Cellpose through the qlivecell interface, producing labeled segmentation masks and associated metadata.

2. **Debris size threshold estimation.**
Estimation of a size-based threshold to distinguish real cells from small segmented objects corresponding to noise or debris.

3. **Removal of segmented debris.**
Cleaning of the segmentation masks by removing objects below the estimated size threshold.

4. **Cell-type counts.**
Quantification of the number of cells in each population within the gastruloid, assuming the presence of multiple cell types.

5. **Local cell density computation.**
Computation of local cell density from the 3D segmentation results, capturing spatial crowding and tissue organization.

6. **Neighborhood composition analysis.**
Extraction of the cellular neighborhood composition for each cell, quantifying the identities of nearby cells. This analysis is only meaningful when more than one cell population is present.

7. **Radial distribution analysis.**
Analysis of the radial distribution of cells relative to the gastruloid center, enabling the study of spatial patterning along the radial axis.

Together, these analyses illustrate how raw 3D image data are transformed into quantitative, single-cell–resolved spatial descriptors of gastruloid organization.


#### Example 2

The second example, found at `examples/example2.ipynb` contains a step by step description of YAP quantification as a nucleus to cytoplasm ratio.

The example dataset consists of a two files in this case:
- `gastruloid2.tif`: A 3D multichannel image of a gastruloid stained for YAP.
- `gastruloid2_background.tif`: A 3D multichannel image of a gastruloid stained for YAP using only the secondary antibody to compute the background signal of such antibody due to inspecific bindings.

The notebook covers the following analysis steps:

1. **Nuclei segmentation and Debris removal.**
As opposed to the first example, in this one segmentation is separated in nucleus vs wholecell and is done using Cellpose. In this first step, the nuclei are segmented and debris is removed using the threshold computed in the previous example.

2. **Whole cell segmentation and Debris removal.**
Segmentation of the whole cells using the membrane label GPI-GFP.

3. **Quantification of YAP background.**
Quantification of YAP in nuclei and cytoplasm in a gastruloid stained only with the YAP secondary antibody. The cytoplasm is computed as the whole-cell masks minus the nuclear mask of that same cell.

4. **Cell classification.**
This dataset lacks the mCherry label, so we used only the emiRFP signal to differentiate the two populations of cells. In this step we quantify a threhold of emiRFP to discern between the two populations.

5. **YAP N/C quantification.**
Quntification of YAP as a Nucleus to Cytoplasm ratio. 

#### Example 3

Finally we present a third example which represent the p53 quantification carried out in the p53 quantification analysis and the p53 degron experiment. 

This example, found at `examples/example3.ipynb` uses a dataset consisting of two files:
- `gastruloid3.tif`: A 3D multichannel image of a gastruloid stained for p53.
- `gastruloid3_background.tif`: A 3D multichannel image of a gastruloid stained for p53 using only the secondary antibody to compute the background signal of such antibody due to inspecific bindings.

The notebook covers the following analysis steps:

1. **Spillover correction.**
Using a global linear calibration model, we correct for signal bleed through from the mCherry and emiRFP channels into the p53 channel.

2. **Segmentation and debris removal.**
Similarly to the first example, we use StarDist for the nuclear segmentation and the same threshold computed in that example for the debris removal

3. **Removal of misclassified cells.**
Detection and removal of cells segmented in the wrong channel. 

4. **p53 quantification.**
Quantification of p53 using the spillover coefficients computed on step one.
---

### Complete Analysis Pipeline (Reference Only)

Within the image-analysis framework, each type of downstream analysis based on the segmentation results is organized into dedicated folders:

- `image-analysis/code/segmentation`
- `image-analysis/code/apoptosis_segmentation`
- `image-analysis/code/counts`
- `image-analysis/code/density`
- `image-analysis/code/neighborhood_composition`
- `image-analysis/code/radial_distribution`
- `image-analysis/code/p53`
- `image-analysis/code/p53_degron`
- `image-analysis/code/YAP`

For additional details, each file contains a short description at the beginning explaining its purpose.