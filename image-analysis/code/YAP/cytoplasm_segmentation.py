"""
Quantify single-cell YAP nuclear-to-cytoplasmic (N/C) ratios.

For each image and condition, loads precomputed whole-cell (GPI-GFP) and nuclear
(DAPI) segmentations, matches nuclei to cells by nearest-centroid tracking, and
derives a cytoplasm mask as (cell mask \ nucleus mask) after overlap QC.

Per matched cell:
- Measure mean nuclear YAP intensity and subtract nuclear background (AB-only).
- Measure mean cytoplasmic YAP intensity and subtract cytoplasmic background.
- Compute N/C = (YAP_nuc_corrected) / (YAP_cyto_corrected).
- Measure nuclear A12 (emiRFP) intensity to classify cells using thresholds
  (low_th / up_th), producing F3-like vs A12-like subsets.

Outputs:
- Per-file CSVs with [Nuclear, Cyto, N/C] for F3-like and A12-like cells.
- Aggregated distributions per condition (WT, KO8, KO25), with outlier removal
  (IQR rule) and KDE + histogram plots saved to figures, plus summary CSVs.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, get_file_names, check_or_create_dir, fill_channels, tif_reader_5D
from qlivecell.celltrack.core.tracking.tracking import greedy_tracking
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import zip_longest
from cellpose import models

model = models.CellposeModel(gpu=True, model_type='cyto3')
channel_names = ["A12", "GPI-GFP", "YAP", "DAPI"]

CONDITIONS = ["WT", "KO8", "KO25", "KO8_ABonly", "KO25_ABonly"]
nuc_quant = []
cyt_quant = []
A12_quant = []

cyt_nuc_ratio = 0.95

YAP_cyt_b = 10255.848867087852
YAP_nuc_b = 10193.577140181687

low_th = 10400
up_th  = 10600

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
}

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/segmentation_results/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/figures/"
path_csvs = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/csvs/"
check_or_create_dir(path_figures)
check_or_create_dir(path_csvs)

for COND in CONDITIONS:
    path_data_dir = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/YAP/2025_02_02_AiryscMultipl_FastMediumQuality_Files/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    
    path_csvs_COND = path_csvs+"{}/".format(COND)
    check_or_create_dir(path_csvs_COND)
    
    ### GET FULL FILE NAME AND FILE CODE ###
    files = get_file_names(path_data_dir)
                    
    ch_cyto = channel_names.index("GPI-GFP")
    ch_nuc = channel_names.index("DAPI")
    ch_yap = channel_names.index("YAP")
    ch_a12 = channel_names.index("A12")

    nuc_quant.append([])
    cyt_quant.append([])
    A12_quant.append([])

    cells_cyto = 0
    cells_nuc = 0

    cells_assigned = 0
    cells_final = 0
    
    for file in files:
        if not ".tif" in file: continue
        
        nuc_quant_file_F3 = []
        cyt_quant_file_F3 = []
        nuc_quant_file_A12 = []
        cyt_quant_file_A12 = []
        
        file, fname = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
        
        path_data = path_data_dir+file
        path_save = path_save_dir+fname
        check_or_create_dir(path_save)

        hyperstack, metadata = tif_reader_5D(path_data)

        chans = fill_channels(ch_nuc, channel_names=channel_names)

        batch_args['name_format'] = "ch"+str(ch_cyto)+"_{}"

        CT_cyto = cellSegTrack(
            path_data,
            path_save,
            batch_args=batch_args,
            channels=chans
        )

        CT_cyto.load()

        batch_args['name_format'] = "ch"+str(ch_nuc)+"_{}"

        CT_nuc = cellSegTrack(
            path_data,
            path_save,
            batch_args=batch_args,
            channels=chans
        )

        CT_nuc.load()
        
        TLabels1 = [cell.label for cell in CT_cyto.jitcells]
        TLabels2 = [cell.label for cell in CT_nuc.jitcells]
        TLabels  = [TLabels1, TLabels2]

        TCenters1 = [cell.centers[0] for cell in CT_cyto.jitcells]
        TCenters2 = [cell.centers[0] for cell in CT_nuc.jitcells]
        TCenters  = [TCenters1, TCenters2]
        
        track_args = {
                "time_step": 1,
                "method": "greedy",
                "dist_th": 7.5,
                "z_th": 2,
            }
        
        xyres = CT_nuc.metadata["XYresolution"]
        zres = CT_nuc.metadata["Zresolution"]

        FinalLabels, label_correspondance = greedy_tracking(TLabels, TCenters, xyres, zres, track_args, lab_max=0)
        
        cells_cyto += len(CT_cyto.jitcells)
        cells_nuc += len(CT_nuc.jitcells)

        cells_assigned = 0
        cells_final = 0
        
        cyt_nuc_ratio = 0.95
        
        for lc in label_correspondance[1]:
            cell1 = CT_cyto._get_cell(lc[1])
            if cell1==None:continue
            cell2 = CT_nuc._get_cell(lc[0])

            z = int(cell1.centers[0][0])
            mask1 = cell1.masks[0][0]
            mask2 = cell2.masks[0][0]
            
            # Assume: mask1 and mask2 are (N, 2) arrays of 2D points
            mask1_view = mask1.view(dtype=[('', mask1.dtype)] * mask1.shape[1])
            mask2_view = mask2.view(dtype=[('', mask2.dtype)] * mask2.shape[1])
            
            cells_assigned += 1
            
            mask_coincidence = np.isin(mask2_view, mask1_view)
            mask_ratio = np.sum(mask_coincidence)/len(mask2)
            if mask_ratio < cyt_nuc_ratio: 
                continue
            
            if len(mask1)/2.0 > len(mask2):
                continue
            
            cells_final +=1
            
            img = hyperstack[0,z,ch_yap]
            nuc_val = np.maximum(np.mean(img[mask2[:, 1], mask2[:, 0]]) - YAP_nuc_b, 0)
            nuc_quant[-1].append(nuc_val)

            img = hyperstack[0,z, ch_a12]
            A12_val = np.mean(img[mask2[:, 1], mask2[:, 0]])
            A12_quant[-1].append(A12_val)

            mask_coincidence_cyto = ~np.isin(mask1_view, mask2_view)
            cyto_mask = np.array([mask1[i] for i, j in enumerate(mask_coincidence_cyto) if j])
            if len(cyto_mask)<5:continue
            img = hyperstack[0,z,ch_yap]
            cyt_val = np.maximum(np.mean(img[cyto_mask[:, 1], cyto_mask[:, 0]]) - YAP_cyt_b, 0)
            cyt_quant[-1].append(cyt_val)

            if A12_val < low_th:
                nuc_quant_file_F3.append(nuc_val)
                cyt_quant_file_F3.append(cyt_val)
            elif A12_val >= up_th:
                nuc_quant_file_A12.append(nuc_val)
                cyt_quant_file_A12.append(cyt_val)

            # from scipy.spatial import ConvexHull
            # r = 200
            # # # cell = np.random.choice(range(len(masks1)))
            # # cell = 23
            # img = hyperstack[0, z, ch_nuc]

            # # mask1 = masks1[cell]
            # # mask2 = masks2[cell]

            # hull = ConvexHull(mask1)
            # outline1 = mask1[hull.vertices]
            # outline1[:] = outline1[:, [1, 0]]

            # hull = ConvexHull(mask2)
            # outline2 = mask2[hull.vertices]
            # outline2[:] = outline2[:, [1, 0]]

            # center = np.mean(outline2, axis=0)

            # fig, ax = plt.subplots()
            # ax.imshow(img)

            # ax.scatter(outline1[:, 1], outline1[:, 0], label="cyto")
            # ax.scatter(outline2[:, 1], outline2[:, 0], label="nuc")

            # ax.set_xlabel("($\mu$m)")
            # ax.set_ylabel("($\mu$m)")

            # xticks = np.arange(0, img.shape[-1]-1, 100)
            # yticks = np.arange(0, img.shape[-1]-1, 100)

            # ax.set_xticks(xticks)
            # ax.set_yticks(yticks)

            # scale = metadata["XYresolution"]
            # ax.set_xticklabels([f"{x*scale:.0f}" for x in xticks])
            # ax.set_yticklabels([f"{y*scale:.0f}" for y in yticks])

            # ax.set_xlim(np.maximum(0, center[1]-r), np.minimum(img.shape[-1], center[1]+r))
            # ax.set_ylim(np.maximum(0, center[0]-r), np.minimum(img.shape[-1], center[0]+r))

            # ax.legend()
            # plt.tight_layout()
            # # plt.savefig("path_figures+cytoplasm.svg", dpi=300)
            # plt.show()

        F3_ratio = np.array(nuc_quant_file_F3)/cyt_quant_file_F3
        A12_ratio = np.array(nuc_quant_file_A12)/cyt_quant_file_A12
        
        # Column names
        col_names = ["Nuclear", "Cyto", "N/C"]

        # Create DataFrame
        df = pd.DataFrame({
            col_names[0]: nuc_quant_file_F3,
            col_names[1]: cyt_quant_file_F3,
            col_names[2]: F3_ratio
        })
        
        path_save_dir_quant = "{}/{}_F3.csv".format(path_csvs_COND, fname)
        
        df.to_csv(path_save_dir_quant, index=False)
        
        # Create DataFrame
        df = pd.DataFrame({
            col_names[0]: nuc_quant_file_A12,
            col_names[1]: cyt_quant_file_A12,
            col_names[2]: A12_ratio
        })
        
        path_save_dir_quant = "{}/{}_A12.csv".format(path_csvs_COND, fname)
        df.to_csv(path_save_dir_quant, index=False)

A12_quant_WT = A12_quant[0]
A12_quant_KO8 = A12_quant[1]
A12_quant_KO25 = A12_quant[2]

cyt_quant_WT = cyt_quant[0]
cyt_quant_KO8 = cyt_quant[1]
cyt_quant_KO25 = cyt_quant[2]

nuc_quant_WT = nuc_quant[0]
nuc_quant_KO8 = nuc_quant[1]
nuc_quant_KO25 = nuc_quant[2]

nuc_quant_WT_WT = [nuc_quant_WT[i] for i in range(len(nuc_quant_WT)) if A12_quant_WT[i]<low_th]
nuc_quant_WT_A12 = [nuc_quant_WT[i] for i in range(len(nuc_quant_WT)) if A12_quant_WT[i]>=up_th]

cyt_quant_WT_WT = [cyt_quant_WT[i] for i in range(len(cyt_quant_WT)) if A12_quant_WT[i]<low_th]
cyt_quant_WT_A12 = [cyt_quant_WT[i] for i in range(len(cyt_quant_WT)) if A12_quant_WT[i]>=up_th]

nuc_quant_KO8_WT = [nuc_quant_KO8[i] for i in range(len(nuc_quant_KO8)) if A12_quant_KO8[i]<low_th]
nuc_quant_KO8_A12 = [nuc_quant_KO8[i] for i in range(len(nuc_quant_KO8)) if A12_quant_KO8[i]>=up_th]

cyt_quant_KO8_WT = [cyt_quant_KO8[i] for i in range(len(cyt_quant_KO8)) if A12_quant_KO8[i]<low_th]
cyt_quant_KO8_A12 = [cyt_quant_KO8[i] for i in range(len(cyt_quant_KO8)) if A12_quant_KO8[i]>=up_th]

nuc_quant_KO25_WT = [nuc_quant_KO25[i] for i in range(len(nuc_quant_KO25)) if A12_quant_KO25[i]<low_th]
nuc_quant_KO25_A12 = [nuc_quant_KO25[i] for i in range(len(nuc_quant_KO25)) if A12_quant_KO25[i]>=up_th]

cyt_quant_KO25_WT = [cyt_quant_KO25[i] for i in range(len(cyt_quant_KO25)) if A12_quant_KO25[i]<low_th]
cyt_quant_KO25_A12 = [cyt_quant_KO25[i] for i in range(len(cyt_quant_KO25)) if A12_quant_KO25[i]>=up_th]

def remove_outliers(data):
    # Calculate the first and third quartiles (Q1 and Q3)
    data_curated = []
    for val in data:
        if not np.isnan(val) and not np.isinf(val):
            data_curated.append(val)
        
    Q1 = np.percentile(data_curated, 25)
    Q3 = np.percentile(data_curated, 75)
    IQR = Q3 - Q1

    # Define the lower and upper bounds for outliers
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    # Filter data to remove outliers
    return [x for x in data_curated if lower_bound <= x <= upper_bound]

from scipy import stats
title = ""
bins=25

xx = np.linspace(0.2, 1.8, 500)

NC_WT = remove_outliers(np.array(nuc_quant_WT_WT)/cyt_quant_WT_WT)
kde_WT = stats.gaussian_kde(NC_WT)

NC_KO8 = remove_outliers(np.array(nuc_quant_KO8_WT)/cyt_quant_KO8_WT)
kde_KO8 = stats.gaussian_kde(NC_KO8)
    
NC_KO25 = remove_outliers(np.array(nuc_quant_KO25_WT)/cyt_quant_KO25_WT)
kde_KO25 = stats.gaussian_kde(NC_KO25)

fig, ax = plt.subplots()
ax.hist(NC_WT, label="F3 N/C WT-WT", alpha=0.3, bins=bins, color="green", density=True)
ax.plot(xx, kde_WT(xx), color="green")

ax.hist(NC_KO8, label="F3 N/C WT-KO8", alpha=0.3, bins=bins, color="cyan", density=True)
ax.plot(xx, kde_KO8(xx), color="cyan")

ax.hist(NC_KO25, label="F3 N/C WT-KO25", alpha=0.3, bins=bins, color="blue", density=True)
ax.plot(xx, kde_KO25(xx), color="blue")

ax.set_xlabel("N/C YAP")
ax.legend()
plt.savefig(path_figures+"F3.svg")
plt.show()

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_WT_WT, cyt_quant_WT_WT, NC_WT, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"WT/F3.csv"
df.to_csv(path_save_dir_quant, index=False)

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_KO8_WT, cyt_quant_KO8_WT, NC_KO8, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"KO8/F3.csv"
df.to_csv(path_save_dir_quant, index=False)

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_KO25_WT, cyt_quant_KO25_WT, NC_KO25, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"KO25/F3.csv"
df.to_csv(path_save_dir_quant, index=False)

NC_WT = remove_outliers(np.array(nuc_quant_WT_A12)/cyt_quant_WT_A12)
kde_WT = stats.gaussian_kde(NC_WT)

NC_KO8 = remove_outliers(np.array(nuc_quant_KO8_A12)/cyt_quant_KO8_A12)
kde_KO8 = stats.gaussian_kde(NC_KO8)

NC_KO25 = remove_outliers(np.array(nuc_quant_KO25_A12)/cyt_quant_KO25_A12)
kde_KO25 = stats.gaussian_kde(NC_KO25)


fig, ax = plt.subplots()
ax.hist(NC_WT, label="A12 N/C WT-WT", alpha=0.3, bins=bins, color="magenta", density=True)
ax.plot(xx, kde_WT(xx), color="magenta")

ax.hist(NC_KO8, label="A12 N/C WT-KO8", alpha=0.3, bins=bins, color="red", density=True)
ax.plot(xx, kde_KO8(xx), color="red")

ax.hist(NC_KO25, label="A12 N/C WT-KO25", alpha=0.3, bins=bins, color="pink", density=True)
ax.plot(xx, kde_KO25(xx), color="pink")

ax.set_xlabel("N/C YAP")
ax.legend()
plt.savefig(path_figures+"A12.svg")
plt.show()

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_WT_A12, cyt_quant_WT_A12, NC_WT, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"WT/A12.csv"
df.to_csv(path_save_dir_quant, index=False)

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_KO8_A12, cyt_quant_KO8_A12, NC_KO8, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"KO8/A12.csv"
df.to_csv(path_save_dir_quant, index=False)

# Create DataFrame
col_names = ["Nuclear", "Cyto", "N/C"]
data = list(zip_longest(nuc_quant_KO25_A12, cyt_quant_KO25_A12, NC_KO25, fillvalue=None))
df = pd.DataFrame(data, columns=col_names)
path_save_dir_quant = path_csvs+"KO25/A12.csv"
df.to_csv(path_save_dir_quant, index=False)


