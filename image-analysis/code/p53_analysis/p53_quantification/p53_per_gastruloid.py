"""
Quantify nuclear p53 per cell and summarize per gastruloid (by z-plane).

For each dataset (WT/KO; n2/n3/n4), loads curated nuclear masks for F3 and A12,
applies spillover correction to the p53 channel (per-marker calibration: slope s
and per-z baseline b0(z)), and computes mean corrected p53 per nucleus. Outputs
per-gastruloid mean p53 per z-plane and the corresponding cell counts for F3 and
A12 (CSV). Also computes WT-weighted z-wise averages (weighted by cell counts)
to normalize F3 p53 profiles and saves the normalized matrix.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=18) 
mpl.rc('axes', labelsize=18) 
mpl.rc('xtick', labelsize=18) 
mpl.rc('ytick', labelsize=18) 
mpl.rc('legend', fontsize=18) 

files_to_exclude = [
    "n2_F3(150)+WT(150)_72h_emiRFP-p53-mCh-DAPI_(40xSil)_Stack1.tif",
    "n2_F3(150)+WT(150)_72h_emiRFP-p53-mCh-DAPI_(40xSil)_Stack2.tif"
]

def build_union_masks(CT_list):
    """
    Build per-z 2D boolean masks marking in-cell pixels from the union of all cells
    across provided CT objects (e.g., CT_F3 and CT_A12).
    Returns: list of length Z with arrays (Y, X) dtype=bool.
    """
    CT0 = CT_list[0]
    Z = CT0.hyperstack.shape[1]
    Y = CT0.hyperstack.shape[-2]
    X = CT0.hyperstack.shape[-1]
    Mz_list = [np.zeros((Y, X), dtype=bool) for _ in range(Z)]
    for CT in CT_list:
        for cell in CT.jitcells:
            z = int(cell.centers[0][0])
            if z < 0 or z >= Z:
                continue
            # find mask for this z
            try:
                zid = cell.zs[0].index(z)
            except ValueError:
                continue
            mask = cell.masks[0][zid]
            yy = mask[:, 1].astype(np.intp)
            xx = mask[:, 0].astype(np.intp)
            Mz_list[z][yy, xx] = True
    return Mz_list

def estimate_b0z_for_file(CT, Mz_list, ch_B, ch_C, s_global, q=0.2):
    # q=0.5 (median) if few high-C cells; q=0.1–0.2 if many might be high
    import numpy as np
    Z = CT.hyperstack.shape[1]
    b0z = np.full(Z, np.nan, dtype=np.float64)
    for z in range(Z):
        Mz = Mz_list[z]
        if not np.any(Mz): continue
        Bz = CT.hyperstack[0, z, ch_B, :, :].astype(np.float64)
        Cz = CT.hyperstack[0, z, ch_C, :, :].astype(np.float64)
        resid = (Cz - s_global * Bz)[Mz].ravel()
        if resid.size < 50: continue
        b0z[z] = float(np.quantile(resid, q))
    # fill empties from available planes
    if np.any(np.isnan(b0z)):
        b0z[np.isnan(b0z)] = np.nanmedian(b0z)
    return b0z


def correct_cell_pixels(CT_ref, mask, z, ch_B, ch_C, s, b0z):
    """Return per-pixel corrected C for one cell at plane z."""
    yy = mask[:, 1].astype(np.intp)
    xx = mask[:, 0].astype(np.intp)
    C_vals = CT_ref.hyperstack[0, z, ch_C, :, :][yy, xx].astype(np.float32)
    B_vals = CT_ref.hyperstack[0, z, ch_B, :, :][yy, xx].astype(np.float32)
    return C_vals - float(b0z[z]) - float(s) * B_vals

CONDS = ["WT", "KO"]
repeats = ["n2", "n3", "n4"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

g_names = []
data_F3 = []
data_A12 = []
cells_F3 = []
cells_A12 = []

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/segmentation_results/"
path_to_save_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/figures/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/p53_per_gastruloid/"
check_or_create_dir(path_to_save_results)

path_to_spillover = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/spillover/"
calibF3 = np.load("{}calibration_F3_to_p53.npz".format(path_to_spillover))
p53_F3_s_global = float(calibF3["s"])
p53_F3_0z = calibF3["b0z"]

calibA12 = np.load("{}calibration_A12_to_p53.npz".format(path_to_spillover))
p53_A12_s_global = float(calibA12["s"])
p53_A12_0z = calibA12["b0z"]

for COND in CONDS:
    for REP in repeats:
        path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/input/{}/{}/".format(COND,REP)
        path_save_dir = "{}{}/{}/".format(master_path_save, COND, REP)
        check_or_create_dir(path_save_dir)
        files = get_file_names(path_data_dir)

        channel_names = ["A12", "p53", "F3", "DAPI"]

        for f, file in enumerate(files):
            
            if file in files_to_exclude: continue

            F3_p53 = [[] for z in range(10)]
            A12_p53 = [[] for z in range(10)]

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
            
            ch_F3 = channel_names.index("F3")
            ch_A12 = channel_names.index("A12")
            ch_p53 = channel_names.index("p53")
            ch_DAPI = channel_names.index("DAPI")

            Mz_list = build_union_masks([CT_F3])
            p53_F3_0z = estimate_b0z_for_file(CT_F3, Mz_list, ch_F3, ch_p53, p53_F3_s_global)
            
            stack_p53 = CT_A12.hyperstack[0,:,ch_p53].astype("float64")
            
            for cell in CT_F3.jitcells:
                center = cell.centers[0]
                z = int(center[0])
                zid = cell.zs[0].index(z)
                mask = cell.masks[0][zid]

                Ccorr_vals = correct_cell_pixels(CT_F3, mask, z, ch_F3, ch_p53, p53_F3_s_global, p53_F3_0z)
                p53_val = float(np.mean(Ccorr_vals))
                F3_p53[z].append(p53_val)
                    
            
            Mz_list = build_union_masks([CT_A12])
            p53_A12_0z = estimate_b0z_for_file(CT_A12, Mz_list, ch_A12, ch_p53, p53_A12_s_global)
            
            for cell in CT_A12.jitcells:
                center = cell.centers[0]
                z = int(center[0])
                zid = cell.zs[0].index(z)
                mask = cell.masks[0][zid]

                Ccorr_vals = correct_cell_pixels(CT_A12, mask, z, ch_A12, ch_p53, p53_A12_s_global, p53_A12_0z)
                p53_val = float(np.mean(Ccorr_vals))
                A12_p53[z].append(p53_val)
            
            
            F3_p53_means = [np.mean(d) for d in F3_p53]
            A12_p53_means = [np.mean(d) for d in A12_p53]
            c_F3 = [len(d) for d in F3_p53]
            c_A12 = [len(d) for d in A12_p53]

            cells_F3.append(c_F3)
            cells_A12.append(c_A12)
            data_F3.append(F3_p53_means)
            data_A12.append(A12_p53_means)
            g_names.append(embcode)
            

import pandas as pd
Z = len(data_F3[0])  # number of z planes, generally 10

df = pd.DataFrame(
    data_F3,
    index=g_names,
    columns=[f"z{j}" for j in range(Z)]
)

# Save to CSV
df.to_csv(path_to_save_results+"F3.csv")

import pandas as pd
Z = len(data_A12[0])  # number of z planes, generally 10

df = pd.DataFrame(
    data_A12,
    index=g_names,
    columns=[f"z{j}" for j in range(Z)]
)

df.to_csv(path_to_save_results+"A12.csv")

import pandas as pd
Z = len(cells_F3[0])  # number of z planes, generally 10

df = pd.DataFrame(
    cells_F3,
    index=g_names,
    columns=[f"z{j}" for j in range(Z)]
)

df.to_csv(path_to_save_results+"F3_number_of_cells.csv")

import pandas as pd
Z = len(cells_A12[0])  # number of z planes, generally 10

df = pd.DataFrame(
    cells_A12,
    index=g_names,
    columns=[f"z{j}" for j in range(Z)]
)

df.to_csv(path_to_save_results+"A12_number_of_cells.csv")

weighted_by_z = []

N_WT = len(g_names[:10])

for z in range(Z):
    means_z = np.array([data_F3[g][z] for g in range(N_WT)])
    weights_z = np.array([cells_F3[g][z] for g in range(N_WT)])
    w_avg_z = np.sum(means_z * weights_z) / np.sum(weights_z)
    weighted_by_z.append(w_avg_z)

N = len(g_names)
data_F3_norm = []
for g in range(N):
    norm_g = np.array(data_F3[g])/weighted_by_z
    data_F3_norm.append(norm_g)
    
import pandas as pd
Z = len(cells_A12[0])  # number of z planes, generally 10

df = pd.DataFrame(
    data_F3_norm,
    index=g_names,
    columns=[f"z{j}" for j in range(Z)]
)
df.to_csv(path_to_save_results+"F3_norm.csv")
