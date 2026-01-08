"""
Estimate cytoplasmic YAP background intensity.

Uses antibody-only control conditions. Loads whole-cell (GPI-GFP) and nuclear (DAPI)
segmentations, matches cells by centroid proximity (greedy tracking), builds a
cytoplasm mask as (whole-cell mask \ nuclear mask), and samples YAP intensity in
that cytoplasm region to obtain the background distribution and threshold.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, get_file_names, tif_reader_5D, check_or_create_dir
import numpy as np
import matplotlib.pyplot as plt

CONDITIONS = ["KO8_ABonly", "KO25_ABonly"]
channel_names = ["A12", "GPI-GFP", "YAP", "DAPI"]

YAP_quant = []

for COND in CONDITIONS:
    path_data_dir = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/YAP/2025_02_02_AiryscMultipl_FastMediumQuality_Files/{}/".format(COND)
    path_save_dir = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/segmentation_results/{}/".format(COND)
    check_or_create_dir(path_save_dir)

    ### GET FULL FILE NAME AND FILE CODE ###
    files = get_file_names(path_data_dir)
                    
    ch_cyto = channel_names.index("GPI-GFP")
    ch_nuc = channel_names.index("DAPI")
    ch_yap = channel_names.index("YAP")
    ch_a12 = channel_names.index("A12")

    nuc_quant = []
    cyt_quant = []
    A12_quant = []

    cells_cyto = 0
    cells_nuc = 0

    cells_assigned = 0
    cells_final = 0

    cyt_nuc_ratio = 0.95
        
    for file in files:
        if not ".tif" in file: continue
        
        file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
        
        path_data = path_data_dir+file
        path_save = path_save_dir+embcode
        check_or_create_dir(path_save)

        hyperstack, metadata = tif_reader_5D(path_data)

        chans = [ch_nuc]
        for _ch in range(len(channel_names)):
            if _ch not in chans:
                chans.append(_ch)

        batch_args = {
            'name_format':"ch"+str(ch_cyto)+"_{}",
            'extension':".tif",
        }

        CT_cyto = cellSegTrack(
            path_data,
            path_save,
            batch_args=batch_args,
            channels=chans
        )
        CT_cyto.load()

        chans = [ch_nuc]
        for _ch in range(len(channel_names)):
            if _ch not in chans:
                chans.append(_ch)

        batch_args = {
            'name_format':"ch"+str(ch_nuc)+"_{}",
            'extension':".tif",
        }

        CT_nuc = cellSegTrack(
            path_data,
            path_save,
            batch_args=batch_args,
            channels=chans
        )
        CT_nuc.load()

        from qlivecell.celltrack.core.tracking.tracking import greedy_tracking
        
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
            
            mask_coincidence = np.isin(mask2_view, mask1_view)
            mask_ratio = np.sum(mask_coincidence)/len(mask2)
            
            if mask_ratio < cyt_nuc_ratio: 
                continue
            
            if len(mask1)/2.0 > len(mask2):
                continue
        
            mask_coincidence_cyto = ~np.isin(mask1_view, mask2_view)
            cyto_mask = np.array([mask1[i] for i, j in enumerate(mask_coincidence_cyto) if j])
            img = hyperstack[0,z,ch_yap]
            YAP_quant.append(np.mean(img[cyto_mask[:, 1], cyto_mask[:, 0]]))

fig, ax = plt.subplots()
ax.hist(YAP_quant, bins=150)
plt.show()

print(np.mean(YAP_quant))
YAP_cy_th = 10255.848867087852
