### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, get_file_names, tif_reader_5D, check_or_create_dir
import numpy as np
import matplotlib.pyplot as plt

CONDITIONS = ["WT", "KO8", "KO25"]
channel_names = ["A12", "GPI-GFP", "YAP", "DAPI"]

A12_quant = []
for COND in CONDITIONS:
    path_data_dir = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/YAP/2025_02_02_AiryscMultipl_FastMediumQuality_Files/{}/".format(COND)
    path_save_dir = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/segmentation_results/{}/".format(COND)
    check_or_create_dir(path_save_dir)

    ### GET FULL FILE NAME AND FILE CODE ###
    files = get_file_names(path_data_dir)
                    
    ch_nuc = channel_names.index("DAPI")
    ch_a12 = channel_names.index("A12")
        
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

        for cell in CT_nuc.jitcells:
            z = int(cell.centers[0][0])
            mask = cell.masks[0][0]
            
            img = hyperstack[0,z, ch_a12]
            A12_quant.append(np.mean(img[mask[:, 1], mask[:, 0]]))

fig, ax = plt.subplots()
ax.hist(A12_quant, bins=150)
plt.show()

# Based on histogram, defined arbitrarily
low_th = 10400
ip_th  = 10600