### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, get_file_names, fill_channels, tif_reader_5D
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy import stats
from qlivecell.celltrack.core.tracking.tracking import greedy_tracking
from cellpose import models
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

model = models.CellposeModel(gpu=True, model_type='cyto3')
channel_names = ["A12", "GPI-GFP", "YAP", "DAPI"]

CONDITIONS = ["WT", "KO8", "KO25"]

cyt_nuc_ratio = 0.95

YAP_cyt_b = 10255.848867087852
YAP_nuc_b = 10193.577140181687

low_th = 10400
up_th  = 10600
                
ch_cyto = channel_names.index("GPI-GFP")
ch_nuc = channel_names.index("DAPI")
ch_yap = channel_names.index("YAP")
ch_a12 = channel_names.index("A12")

cells_cyto = 0
cells_nuc = 0

cells_assigned = 0
cells_final = 0

COND = CONDITIONS[0]

path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/figures/"
path_data_dir = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/YAP/2025_02_02_AiryscMultipl_FastMediumQuality_Files/{}/".format(COND)
path_save_dir = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/segmentation_results/{}/".format(COND)

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
}

files = get_file_names(path_data_dir) 
files = [file for file in files if ".tif" in file]   
file = files[2]    
file, fname = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)

path_data = path_data_dir+file
path_save = path_save_dir+fname

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

chans = fill_channels(ch_nuc, channel_names=channel_names)
batch_args['name_format'] = "ch"+str(ch_nuc)+"_{}"
CT_nuc = cellSegTrack(
    path_data,
    path_save,
    batch_args=batch_args,
    channels=chans
)
CT_nuc.load()

CT_nuc.plot_tracking()

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

nuc_quant_file_F3 = []
cyt_quant_file_F3 = []
nuc_quant_file_A12 = []
cyt_quant_file_A12 = []

masks1 = []
masks2 = []
masks3 = []
zs = []
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


    img = hyperstack[0,z, ch_a12]
    A12_val = np.mean(img[mask2[:, 1], mask2[:, 0]])

    mask_coincidence_cyto = ~np.isin(mask1_view, mask2_view)
    cyto_mask = np.array([mask1[i] for i, j in enumerate(mask_coincidence_cyto) if j])
    img = hyperstack[0,z,ch_yap]
    cyt_val = np.maximum(np.mean(img[cyto_mask[:, 1], cyto_mask[:, 0]]) - YAP_cyt_b, 0)

    if A12_val < low_th:
        nuc_quant_file_F3.append(nuc_val)
        cyt_quant_file_F3.append(cyt_val)
    elif A12_val >= up_th:
        nuc_quant_file_A12.append(nuc_val)
        cyt_quant_file_A12.append(cyt_val)

    masks1.append(mask1)
    masks2.append(mask2)
    masks3.append(cyto_mask)
    zs.append(z)

F3_ratio = np.array(nuc_quant_file_F3)/cyt_quant_file_F3
A12_ratio = np.array(nuc_quant_file_A12)/cyt_quant_file_A12

def remove_outliers(data):
    # Calculate the first and third quartiles (Q1 and Q3)
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1

    # Define the lower and upper bounds for outliers
    lower_bound = Q1 - 2.5 * IQR
    upper_bound = Q3 + 2.5 * IQR

    # Filter data to remove outliers
    return [x for x in data if lower_bound <= x <= upper_bound]

title = ""
bins=25

xx = np.linspace(0.2, 1.8, 500)

F3 = []
for val in F3_ratio:
    if not np.isnan(val) and not np.isinf(val):
        F3.append(val)

A12 = []
for val in A12_ratio:
    if not np.isnan(val) and not np.isinf(val):
        A12.append(val)
        
F3 = remove_outliers(F3)
kde_F3 = stats.gaussian_kde(F3)

A12 = remove_outliers(A12)
kde_A12 = stats.gaussian_kde(A12)

fig, ax = plt.subplots()
ax.hist(F3, label="F3 N/C", alpha=0.3, bins=bins, color="green", density=True)
ax.plot(xx, kde_F3(xx), color="green")

ax.hist(A12, label="A12 N/C", alpha=0.3, bins=bins, color="magenta", density=True)
ax.plot(xx, kde_A12(xx), color="magenta")

ax.set_xlabel("N/C YAP")
ax.legend()
plt.savefig(path_figures+"cyto_quantification.svg")
plt.show()

r = 200
# cell = np.random.choice(range(len(masks1)))
cell = 23
img = hyperstack[0, zs[cell], ch_nuc]

mask1 = masks1[cell]
mask2 = masks2[cell]

hull = ConvexHull(mask1)
outline1 = mask1[hull.vertices]
outline1[:] = outline1[:, [1, 0]]

hull = ConvexHull(mask2)
outline2 = mask2[hull.vertices]
outline2[:] = outline2[:, [1, 0]]

center = np.mean(outline2, axis=0)

fig, ax = plt.subplots()
ax.imshow(img)

ax.scatter(outline1[:, 1], outline1[:, 0], label="cyto")
ax.scatter(outline2[:, 1], outline2[:, 0], label="nuc")

ax.set_xlabel("($\mu$m)")
ax.set_ylabel("($\mu$m)")

xticks = np.arange(0, img.shape[-1]-1, 100)
yticks = np.arange(0, img.shape[-1]-1, 100)

ax.set_xticks(xticks)
ax.set_yticks(yticks)

scale = metadata["XYresolution"]
ax.set_xticklabels([f"{x*scale:.0f}" for x in xticks])
ax.set_yticklabels([f"{y*scale:.0f}" for y in yticks])

ax.set_xlim(np.maximum(0, center[1]-r), np.minimum(img.shape[-1], center[1]+r))
ax.set_ylim(np.maximum(0, center[0]-r), np.minimum(img.shape[-1], center[0]+r))

ax.legend()
plt.tight_layout()
plt.savefig("path_figures+cytoplasm.svg", dpi=300)
plt.show()
