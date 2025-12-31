# In this file we compute the debris threshold for each gastruloid

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=14) 
mpl.rc('axes', labelsize=14) 
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('legend', fontsize=14) 

from scipy.signal import argrelextrema
from sklearn.neighbors import KernelDensity

# ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
CONDS = ["auxin_48-72_48", "auxin_48-72_72" , "auxin_48-72_96", "auxin_72-96_72", "auxin_72-96_96", "noauxin_72", "noauxin_96", "secondaryonly"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

thresholds = []

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53_degron/segmentation_results/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53_degron/figures/"
check_or_create_dir(path_to_save_results)

areas = []
for COND in CONDS:
    path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/2025_09_09_OsTIRMosaic_p53Timecourse/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    
    files = get_file_names(path_data_dir)
    
    channel_names = ["A12", "p53", "F3", "DAPI"]

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
        
        for cell in CT_F3.jitcells:
            zc = int(cell.centers[0][0])
            zcid = cell.zs[0].index(zc)

            mask = cell.masks[0][zcid]
            area = len(mask)* CT_F3.CT_info.xyresolution**2
            areas.append(area)

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

        for cell in CT_A12.jitcells:
            zc = int(cell.centers[0][0])
            zcid = cell.zs[0].index(zc)

            mask = cell.masks[0][zcid]
            area = len(mask)* CT_A12.CT_info.xyresolution**2
            areas.append(area)


fig, ax = plt.subplots(2,1, figsize=(5,8))
threshold = 0
data = np.array(areas)
ax[0].hist(data, bins=200, color=[0.0, 0.8, 0.0], density=True, alpha=0.6,)

x = np.arange(0, step=0.1, stop=np.max(data))
bw = 7
modelo_kde = KernelDensity(kernel="linear", bandwidth=bw)
modelo_kde.fit(X=data.reshape(-1, 1))
densidad_pred = np.exp(modelo_kde.score_samples(x.reshape((-1, 1))))
ax[0].plot(x, densidad_pred, color="magenta")

local_minima = argrelextrema(densidad_pred, np.less)[0]
threshold = x[local_minima[0]]
x_th = np.ones(len(x)) * x[local_minima[0]]
y_th = np.linspace(0, np.max(densidad_pred), num=len(x))
ax[0].plot(x_th, y_th, c="k", ls="--",lw=2, label="debris th.")

thresholds.append(threshold)

ax[0].set_ylabel("count")
ax[1].set_ylabel("count")

ax[1].hist(data, bins=200, color=[0.0, 0.8, 0.0], density=True, alpha=0.6)
ax[1].set_xlabel(r"area ($\mu$m$^2$)")
ax[1].plot(x, densidad_pred, color="magenta")
ax[1].plot(x_th, y_th, c="k", ls="--",lw=2, label="debris th.")

ax[1].set_xlim(-1, 75)
plt.tight_layout()

plt.savefig(path_to_save_results+"debris_thresholds.pdf")
plt.show()
