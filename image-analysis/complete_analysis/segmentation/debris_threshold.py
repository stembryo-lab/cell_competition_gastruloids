"""
Compute debris-removal size thresholds from segmentation area distributions.

Loads existing 3D segmentations, collects per-object center areas (µm²) across
embryos/conditions, fits a KDE to the area distribution, and uses the first local
minimum of the density as the debris threshold. Outputs summary plots to PDF.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
CONDS = ["WT", "KO"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
master_path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/debris_removal/"

areas = {}
for EXP in EXPERIMENTS:
    if EXP=="2023_11_17_Casp3":
        TIMES = ["48hr", "72hr", "96hr"]
    elif EXP=="2024_03_Casp3":        
        TIMES = ["48hr", "60hr", "72hr", "96hr"]
    
    areas[EXP] = {}

    for TIME in TIMES:
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if TIME=="96hr":
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        elif EXP=="2024_03_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
    
        areas[EXP][TIME] = []
        for COND in CONDS:
            path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
            path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
            files = get_file_names(path_data_dir)
            

            for f, file in enumerate(files):
                path_data = path_data_dir+file
                file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
                path_save = path_save_dir+embcode
                check_or_create_dir(path_save)
                
                ch = channel_names.index("F3")
                
                batch_args['name_format'] = "ch"+str(ch)+"_{}"
                
                chans = fill_channels(channel=ch, channel_names=channel_names)
                        
                CT_F3 = cellSegTrack(
                    path_data,
                    path_save,
                    batch_args=batch_args,
                    channels=chans
                )

                CT_F3.load()
                    
                ch = channel_names.index("A12")
                batch_args['name_format'] = "ch"+str(ch)+"_{}"                

                chans = fill_channels(channel=ch, channel_names=channel_names)

                CT_A12 = cellSegTrack(
                    path_data,
                    path_save,
                    batch_args=batch_args,
                    channels=chans
                )

                CT_A12.load()

                for cell in CT_F3.jitcells:
                    zc = int(cell.centers[0][0])
                    zcid = cell.zs[0].index(zc)

                    mask = cell.masks[0][zcid]
                    area = len(mask) * CT_F3.metadata["XYresolution"]**2
                    areas[EXP][TIME].append(area)
                    
                for cell in CT_A12.jitcells:
                    zc = int(cell.centers[0][0])
                    zcid = cell.zs[0].index(zc)

                    mask = cell.masks[0][zcid]
                    area = len(mask) * CT_A12.metadata["XYresolution"]**2
                    areas[EXP][TIME].append(area)
                        
TIMES = ["48hr", "72hr", "96hr"]
fig, ax = plt.subplots(2, 3,figsize=(12,6))
for e, EXP in enumerate(EXPERIMENTS):
    for t, TIME in enumerate(TIMES):
        if t==0:
            ax[e, t].set_ylabel(EXP)
        if e==0:
            ax[e, t].set_title(TIME)
        threshold = 0   
        data = np.array(areas[EXP][TIME])
        data = np.clip(data, 0, 250)
        x = np.arange(0, step=0.5, stop=np.max(data))
        bw = 15
        modelo_kde = KernelDensity(kernel="linear", bandwidth=bw)
        modelo_kde.fit(X=data.reshape(-1, 1))
        densidad_pred = np.exp(modelo_kde.score_samples(x.reshape((-1, 1))))
        local_minima = argrelextrema(densidad_pred, np.less)[0]
        threshold = x[local_minima[0]]
        x_th = np.ones(len(x)) * x[local_minima[0]]
        y_th = np.linspace(0, np.max(densidad_pred), num=len(x))
        
        ax[e, t].hist(data, bins=100, color=[0.0, 0.8, 0.0], density=True, alpha=0.6)
        ax[e, t].plot(x, densidad_pred, color="magenta")
        ax[e, t].plot(x_th, y_th, c="k", ls="--",lw=2, label="debris th. = {:0.1f}".format(threshold))
        ax[e, t].set_xlim(-1,200)
        ax[e, t].set_xlabel("area (µm²)")
        ax[e, t].set_yticks([])
        ax[e, t].legend(loc="upper left")
plt.tight_layout()
plt.savefig(master_path_to_save_results+"debris_thresholds.pdf")
plt.show()

a = "{:0.1f}".format(threshold)