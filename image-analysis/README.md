# Image analysis

Image analysis for this work was primarily focused on using cellular segmentation to extract spatial information in gastruloids at the single-cell level.  
For this purpose, we used the software **qlivecell**, available at:  
https://github.com/dsb-lab/qlivecell

In particular, we used version `v0.8`. To install this version, run:

```bash
python -m pip install 'qlivecell @ git+https://github.com/dsb-lab/qlivecell@0.8'
```

The `qlivecell` package provides tools for segmentation, error correction, and downstream analysis.  
Cellular segmentation relies on **StarDist** (https://github.com/stardist/stardist) and **Cellpose** (https://github.com/MouseLand/cellpose). Please refer to their respective GitHub pages for installation instructions.

Inside the `cell_competition_gastruloids/image-analysis` folder you will find two subfolders:

- `examples`  
- `complete_analysis`

---

#### `examples` folder

The `examples` folder contains minimal, runnable examples of all analyses performed in this work, each using a single gastruloid.  
These scripts are designed to be executable and to illustrate the analysis workflow step by step.

#### `complete_analysis` folder

The `complete_analysis` folder contains all the code used to generate the image-analysis results presented in the article.  
This folder is provided **for reference only**, as the scripts will not run directly due to the absence of the full dataset in this repository.

Below, you will find:

1. A complete guide on how to run the example scripts.  
2. A walkthrough of the structure of the `complete_analysis` folder, allowing the reader to map each script to its corresponding figure or panel in the paper.

---

### Example Workflows (Runnable)

The `examples` folder provides a minimal, end-to-end workflow illustrating how 3D gastruloid images are processed, segmented, and analyzed using **qlivecell**.  
All examples operate on a single gastruloid image and are intended to be run sequentially, as later steps depend on the outputs of earlier ones.

The example dataset consists of a single file:

- `gastruloid.tif`: A 3D multichannel image of a gastruloid used as a proof of concept for the full analysis pipeline.

Before running the analyses, make sure that `qlivecell` is installed.  
Each script can be executed using:

```bash
python script_name.py
```

Before running the scripts, ensure that the path to the example gastruloid is correctly specified.  
Below, we describe the role of each script. Scripts are listed in the recommended execution order.

---

#### 1. Cell segmentation

**Script:** `segmentation.py`

This script performs the initial 3D cell segmentation of the gastruloid using **StarDist** and **Cellpose** through the `qlivecell` interface.  
It produces labeled segmentation masks and associated metadata that are used as input for all downstream analyses.

Typical outputs include:
- 3D cell masks  
- Cell centroids  
- Per-cell geometric properties  

This step is mandatory for all subsequent scripts.

---

#### 2. Debris size threshold estimation

**Script:** `debris_threshold.py`

Segmentation may include small objects that do not correspond to real cells (e.g. noise or debris).  
This script illustrates how to estimate a size-based threshold to distinguish cells from debris based on the distribution of segmented object sizes.

The computed threshold is used in the next step to clean the segmentation.

---

#### 3. Removal of segmented debris

**Script:** `debris_removal.py`

Using the threshold obtained in the previous step, this script removes small segmented objects that are unlikely to represent real cells.  
The resulting cleaned segmentation masks are used for all quantitative analyses.

---

#### 4. Cell-type counts

**Script:** `counts.py`

This example shows how to use the cleaned segmentation results to count cells from each population within the gastruloid.  
It assumes the presence of multiple cell types (e.g. different fluorescent nuclear markers) and illustrates how to extract population-specific counts.

---

#### 5. Nuclear fluorescence quantification

**Script:** `nuclear_quantification.py`

This script demonstrates how to quantify nuclear fluorescence at the single-cell level using the segmentation masks.  
It illustrates how to:
- Use nuclear masks  
- Perform background subtraction  
- Extract per-cell fluorescence intensities  

The resulting values can be used for downstream statistical or spatial analyses.

---

#### 6. Local cell density computation

**Script:** `density.py`

This example computes local cell density from the 3D segmentation results.  
It illustrates how spatial relationships between cells are used to estimate density metrics, which are later employed to study tissue organization and crowding effects.

---

#### 7. Neighborhood composition analysis

**Script:** `neighborhood_composition.py`

This script extracts the cellular neighborhood composition for each cell in the gastruloid.  
For every cell, the identities of nearby cells are quantified, enabling analysis of local cellular environments.

**Note:** This analysis is only meaningful when more than one cell population is present in the gastruloid.

---

#### 8. Radial distribution analysis

**Script:** `radial_distribution.py`

Finally, this example computes the radial distribution of cells relative to the gastruloid center.  
It shows how spatial coordinates derived from the segmentation are used to analyze spatial patterning along the radial axis of the gastruloid.

---

Together, these scripts illustrate how raw 3D image data are transformed into quantitative, single-cell–resolved spatial descriptors of gastruloid organization.

---

### Full Analysis Pipeline (Reference Only)

Within the image-analysis framework, each type of downstream analysis based on the segmentation results is organized into dedicated folders:

- `image-analysis/code/segmentation`: Code to perform segmentation and post-processing of segmentation results  
- `image-analysis/code/apoptosis_segmentation`: Code for segmentation of apoptotic events using Caspase-3 stainings  
- `image-analysis/code/counts`: Code to generate cell and apoptotic counts  
- `image-analysis/code/density`: Code for density analyses  
- `image-analysis/code/neighborhood_composition`: Code for neighborhood composition analyses  
- `image-analysis/code/radial_distribution`: Code for radial distribution analyses  
- `image-analysis/code/p53`: Code for p53 quantification  
- `image-analysis/code/p53_degron`: Code for the p53 degron experiment  
- `image-analysis/code/YAP`: Code for YAP quantification  

For additional details, each file contains a short description at the beginning explaining its purpose.