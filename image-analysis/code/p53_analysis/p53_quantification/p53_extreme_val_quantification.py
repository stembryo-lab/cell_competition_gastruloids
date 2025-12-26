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

# 4.5 iqr
extreme_val_thresholds = [8136.314674377441, 6006.656898498535, 5399.654033660889, 4638.458486557007, 4094.4220390319824, 3730.2568321228027, 3143.496006011963, 2677.543336868286, 2354.30233001709, 1920.6456594467163]

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/segmentation_results/"
path_to_save_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/figures/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/p53_extreme_quantification/"
check_or_create_dir(path_to_save_results)

path_to_spillover = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/spillover/"
calibF3 = np.load("{}calibration_F3_to_p53.npz").format(path_to_spillover)
p53_F3_s_global = float(calibF3["s"])
p53_F3_0z = calibF3["b0z"]

calibA12 = np.load("{}calibration_A12_to_p53.npz").format(path_to_spillover)
p53_A12_s_global = float(calibA12["s"])
p53_A12_0z = calibA12["b0z"]

ExtremesF3 = {}
ExtremesA12 = {}
for COND in CONDS:
    ExtremesF3[COND] = {}
    ExtremesA12[COND] = {}

    for REP in repeats:
        ExtremesF3[COND][REP] = []
        ExtremesA12[COND][REP] = []

        path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/input/{}/{}/".format(COND,REP)
        path_save_dir = "{}{}/{}/".format(master_path_save, COND, REP)
        check_or_create_dir(path_save_dir)
        files = get_file_names(path_data_dir)

        channel_names = ["A12", "p53", "F3", "DAPI"]
                    
        ch_F3 = channel_names.index("F3")
        ch_A12 = channel_names.index("A12")
        ch_p53 = channel_names.index("p53")
        ch_DAPI = channel_names.index("DAPI")
        
        for f, file in enumerate(files):
            
            if file in files_to_exclude: continue
            ExtremesF3[COND][REP].append(0)
            ExtremesA12[COND][REP].append(0)
            
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
            
            zn = CT_F3.hyperstack.shape[1]
            p53_F3 = [[] for z in range(zn)]
            p53_A12 = [[] for z in range(zn)]

            for cell in CT_F3.jitcells:
                center = cell.centers[0]
                z = int(center[0])
                zid = cell.zs[0].index(z)
                mask = cell.masks[0][zid]

                Ccorr_vals = correct_cell_pixels(CT_F3, mask, z, ch_F3, ch_p53, p53_F3_s_global, p53_F3_0z)
                p53_val = float(np.mean(Ccorr_vals))
                p53_F3[z].append(p53_val)
                if p53_val > extreme_val_thresholds[z]:
                    ExtremesF3[COND][REP][-1]+=1
            
            ExtremesF3[COND][REP][-1]/=len(CT_F3.jitcells)
            
            for cell in CT_A12.jitcells:
                center = cell.centers[0]
                z = int(center[0])
                zid = cell.zs[0].index(z)
                mask = cell.masks[0][zid]

                Ccorr_vals = correct_cell_pixels(CT_A12, mask, z, ch_A12, ch_p53, p53_A12_s_global, p53_A12_0z)
                p53_val = float(np.mean(Ccorr_vals))
                p53_A12[z].append(p53_val)
                if p53_val > extreme_val_thresholds[z]:
                    ExtremesA12[COND][REP][-1]+=1
                    
            ExtremesA12[COND][REP][-1]/=len(CT_A12.jitcells)
            

ExtremesF3_WT_means = [np.mean(ExtremesF3["WT"][REP]) for REP in repeats]
ExtremesF3_WT_stds = [np.std(ExtremesF3["WT"][REP]) for REP in repeats]

ExtremesF3_KO_means = [np.mean(ExtremesF3["KO"][REP]) for REP in repeats]
ExtremesF3_KO_stds = [np.std(ExtremesF3["KO"][REP]) for REP in repeats]

# Set up bar positions
x = np.arange(len(repeats))   # one position per REP
width = 0.35                  # width of each bar

fig, ax = plt.subplots()

# Bars
bars_WT = ax.bar(x - width/2, ExtremesF3_WT_means, width,
                 yerr=ExtremesF3_WT_stds, label="F3 with WT",
                 color=(0.0, 0.8, 0.0), capsize=5)

bars_KO = ax.bar(x + width/2, ExtremesF3_KO_means, width,
                 yerr=ExtremesF3_KO_stds, label="F3 with KO",
                 color="cyan", capsize=5)

# Labels in the middle of each pair
ax.set_xticks(x)                  # one tick per pair
ax.set_xticklabels(repeats)       # label with REP names

# Axes labels & legend
ax.set_ylabel("Proportion of p53-high cells")
ax.set_title("WT vs KO (F3)")
ax.legend()

plt.show()

all_extremeF3_WT = [val for REP in repeats for val in ExtremesF3["WT"][REP]]
all_extremeF3_KO = [val for REP in repeats for val in ExtremesF3["KO"][REP]]

all_extremeA12_WT = [val for REP in repeats for val in ExtremesA12["WT"][REP]]
all_extremeA12_KO = [val for REP in repeats for val in ExtremesA12["KO"][REP]]


from scipy.stats import ttest_ind
from scipy.stats import shapiro

stat, p_shapiro = shapiro(all_extremeF3_WT)
print("WT Shapiro-Wilk p =", p_shapiro)

stat, p_shapiro = shapiro(all_extremeF3_KO)
print("KO Shapiro-Wilk p =", p_shapiro)

t_stat, p_t = ttest_ind(all_extremeF3_WT, all_extremeF3_KO, equal_var=False)
print(f"T-test: t={t_stat:.3f}, p={p_t:.3e}")

from scipy.stats import mannwhitneyu

u_stat, p_mw = mannwhitneyu(all_extremeF3_WT, all_extremeF3_KO, alternative="two-sided")
print(f"Mann–Whitney U: U={u_stat}, p={p_mw:.3e}")

from scipy.stats import ks_2samp

ks_stat, p_ks = ks_2samp(all_extremeF3_WT, all_extremeF3_KO)
print(f"KS test: D={ks_stat:.3f}, p={p_ks:.3e}")

# Bar positions
labels = ["A12 - KO","A12 - WT", "F3 with WT", "F3 with KO"]
means = [np.mean(all_extremeA12_KO), np.mean(all_extremeA12_WT), np.mean(all_extremeF3_WT), np.mean(all_extremeF3_KO)]
stds = [np.std(all_extremeA12_KO), np.std(all_extremeA12_WT), np.std(all_extremeF3_WT), np.std(all_extremeF3_KO)]
colors = [(0.6, 0.0, 0.6), (0.8, 0.0, 0.8), (0.0, 0.8, 0.0), "cyan"]

x = np.arange(len(labels))

# Plot
fig, ax = plt.subplots()
bars = ax.bar(x, means, yerr=stds, capsize=5, color=colors)

# Labels and formatting
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylabel("Proportion of p53-high cells")
plt.show()

bins = 10
fig, ax = plt.subplots()
ax.hist(all_extremeA12_KO, bins=bins, color=colors[0], alpha=0.5)
ax.hist(all_extremeA12_WT, bins=bins, color=colors[1], alpha=0.5)
ax.hist(all_extremeF3_WT, bins=bins, color=colors[2], alpha=0.5)
ax.hist(all_extremeF3_KO, bins=bins, color=colors[3], alpha=0.5)
plt.show()

from scipy.stats import ttest_ind

# Define the comparisons you want to annotate (pairs of bar indices)
comparisons = [(1,2), (2,3)]  # e.g. A12 KO vs A12 WT, F3 WT vs F3 KO

# Bar positions
labels = ["A12 - KO", "A12 - WT", "F3 with WT", "F3 with KO"]
data = [all_extremeA12_KO, all_extremeA12_WT, all_extremeF3_WT, all_extremeF3_KO]

means = [np.mean(d) for d in data]
stds  = [np.std(d) for d in data]
colors = [(0.6, 0.0, 0.6), (0.8, 0.0, 0.8), (0.0, 0.8, 0.0), "cyan"]
x = np.arange(len(labels))

# Run tests and add bars
ymax = max(max(d) for d in data) * 1.01  # starting height for bars
h = 0.02 * ymax   # height of bars
step = 0

fig, ax = plt.subplots(figsize=(8,5))

# Plot bars
bars = ax.bar(x, means, yerr=stds, capsize=5, color=colors, alpha=0.7)

# Overlay scatter points
for i, vals in enumerate(data):
    # random jitter inside the bar width
    jitter = np.random.uniform(-0.2, 0.2, size=len(vals))  
    ax.scatter(np.full(len(vals), x[i]) + jitter, vals,
               color="black", s=20, alpha=0.7, zorder=3)

# Labels and formatting
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=20)
ax.set_ylabel("Proportion of p53-high cells")
ax.set_title("Bar plot with individual datapoints")

for i, j in comparisons:
    vals1, vals2 = data[i], data[j]

    # Example: Welch’s t-test
    _, p = ttest_ind(vals1, vals2, equal_var=False)

    # Bar coordinates
    x1, x2 = x[i], x[j]
    y = ymax + step*h
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='k')

    # Decide significance stars
    if p < 0.001:
        stars = '***'
    elif p < 0.01:
        stars = '**'
    elif p < 0.05:
        stars = '*'
    else:
        stars = 'n.s.'  # not significant

    ax.text((x1+x2)/2, y+h, stars, ha='center', va='bottom', fontsize=12)

    step += 4  # space out multiple bars

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig("{}barplots.svg".format(path_to_save_figures))
plt.show()

import pandas as pd

# Keep labels and data in the same order you plot
labels = ["A12 - KO", "A12 - WT", "F3 with WT", "F3 with KO"]
data = [all_extremeA12_KO, all_extremeA12_WT, all_extremeF3_WT, all_extremeF3_KO]

# Build a dict of Series to handle different lengths (NaN padding)
cols = {lab: pd.Series(vals) for lab, vals in zip(labels, data)}

df = pd.DataFrame(cols)
# Save to CSV
df.to_csv(path_to_save_results+"barplot_underlying_data.csv", index=False)

