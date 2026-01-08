"""
Visual example of the radial-position metric on a single gastruloid.

Loads saved segmentations, segments the gastruloid edge on one z-slice, and
illustrates the metric by plotting a representative cell, the gastruloid centroid,
and the closest point on the edge (with the corresponding distance vectors).
Exports the figure as SVG/PDF.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels, compute_dists_jit, compute_distance_xy_jit, EmbryoSegmentation, tif_reader_5D
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

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
path_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/values/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/figures/"

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
EXP = EXPERIMENTS[0]
path_results_exp = path_results+"{}/".format(EXP)
path_figures_exp = path_figures+"{}/".format(EXP)
check_or_create_dir(path_figures_exp)

TIMES = ["48hr", "72hr", "96hr"]
TIME = TIMES[0]

if EXP=="2023_11_17_Casp3":
    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
    if TIME=="96hr":
        channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
elif EXP=="2024_03_Casp3":
    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
    
CONDS = ["WT", "KO"]
COND = CONDS[1]

apo_stages = ["early", "mid", "late"]
apo_stage = apo_stages[0]
 
path_save_results='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, COND)
path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
files_data = get_file_names(path_data_dir)
files = get_file_names(path_data_dir)
file = files[0]

path_data = path_data_dir+file
file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
path_save = correct_path(path_save_dir+embcode)

ch_F3 = channel_names.index("F3")
ch_A12 = channel_names.index("A12")
ch_Casp3 = channel_names.index("Casp3")
            
batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

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

ch = channel_names.index("Casp3")
batch_args['name_format'] = "ch"+str(ch)+"_{}_"+apo_stage                
chans = fill_channels(channel=ch, channel_names=channel_names)

CT_Casp3 = cellSegTrack(
    path_data,
    path_save,
    batch_args=batch_args,
    channels=chans
)

CT_Casp3.load()

points = []
zs = []
for cell in CT_F3.jitcells:
    outlines = cell.outlines[0]
    for zid, z in enumerate(cell.zs[0]):
        for point in outlines[zid]:
            point3D = [z, point[0], point[1]]
            points.append(point3D)

for cell in CT_A12.jitcells:
    outlines = cell.outlines[0]
    for zid, z in enumerate(cell.zs[0]):
        for point in outlines[zid]:
            point3D = [z, point[0], point[1]]
            points.append(point3D)


xyres = CT_F3.metadata["XYresolution"]
zres  = CT_F3.metadata["Zresolution"]
resolutions = np.array([zres, xyres, xyres])
points = np.array(points)*resolutions

from scipy.spatial import ConvexHull
hull = ConvexHull(points)
outline3D = points[hull.vertices]

# Get embryo 3D Centroid
centers = []
for cell in CT_F3.jitcells:
    center = cell.centers[0]
    centers.append(center)
    
for cell in CT_A12.jitcells:
    center = cell.centers[0]
    centers.append(center)

centers = np.array(centers)
minz = int(np.min(centers[:,0]))
maxz = int(np.max(centers[:,0]))

centers = np.array(centers)*resolutions
centroid = np.mean(centers, axis=0)

from qlivecell import EmbryoSegmentation, tif_reader_5D
hyperstack, metadata = tif_reader_5D(path_data)
channels_seg = np.array([ch_A12, ch_F3, ch_Casp3])
hyperstack_seg = np.sum(hyperstack[:,:,channels_seg, :, :].astype("int32"), axis=2)

z_plot = np.rint(hyperstack_seg.shape[1]/2).astype("int64")
ES = EmbryoSegmentation(
        hyperstack_seg,
        ksize=5,
        ksigma=20,
        binths=8,
        apply_biths_to_zrange_only=False,
        checkerboard_size=10,
        num_inter=100,
        smoothing=20,
        trange=None,
        zrange=range(minz, maxz+1),
        mp_threads=14,
    )

ES(hyperstack_seg)

plt.imshow(ES.LS[0][z_plot])
plt.show()

import numpy as np
from skimage import measure

contour_points3D = []
for zid, z in enumerate(range(z_plot, z_plot+1)):
    print(z)
    contours = measure.find_contours(ES.LS[0][z], 0.5)
    contour = []
    for cont in contours:
        if len(cont)>len(contour):
            contour = cont
    for p in contour:
        contour_points3D.append(np.array([p[1], p[0]]))
contour_points3D = np.array(contour_points3D)

print("got contours")

centers_all = []
for cell in CT_F3.jitcells:
    if z_plot in cell.zs[0]:
        zid = cell.zs[0].index(z_plot)
        centers_all.append(cell.centers_all[0][zid][1:])

for cell in CT_A12.jitcells:
    if z_plot in cell.zs[0]:
        zid = cell.zs[0].index(z_plot)
        centers_all.append(cell.centers_all[0][zid][1:])
centers_all = np.array(centers_all)

point_inside=centers_all[2]
dists = compute_dists_jit(np.array([point_inside]), np.array(contour_points3D), compute_distance_xy_jit)
closest = np.argmin(dists, axis=1)
point_contour = contour_points3D[closest][0]
centroid = np.mean(centers_all, axis=0)

fig, ax = plt.subplots(figsize=(7,7))
# Plot the ellipse
ax.imshow(hyperstack_seg[0, z_plot])
ax.plot(contour_points3D[:, 0], contour_points3D[:, 1], ls='-', label='edge', c='k', lw=4)
ax.plot([point_inside[0], point_contour[0]], [point_inside[1], point_contour[1]], c=[0.8,0.8,0.8], lw=3)
ax.plot([point_inside[0], centroid[0]], [point_inside[1], centroid[1]], c=[0.3,0.3,0.3], lw=3)
ax.scatter([point_inside[0]], [point_inside[1]], s=75, label="point", color="brown", zorder=10)
ax.scatter([centroid[0]], [centroid[1]], s=75, label="centroid", color="k", zorder=10)
ax.scatter([point_contour[0]], [point_contour[1]], s=75, label="closest edge", color="yellow", zorder=10)

ax.set_title('relative position on gastruloid')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1.0))
ax.set_aspect('equal')
ax.spines[['bottom','left', 'right', 'top']].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(100, 512-100)
ax.set_ylim(100, 512-100)
plt.tight_layout()
plt.savefig(path_figures_exp+"radial_example.svg")
plt.savefig(path_figures_exp+"radial_example.pdf")
plt.show()

