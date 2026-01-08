"""
Compare the radial distribution of apoptotic cells to the full cell population.

Loads precomputed radial-distance arrays and computes the relative radial position
d = d_centroid / (d_centroid + d_edge) for all cells vs. Casp3+ cells (early/mid/late).
Generates summary plots across time points (WT and KO) and saves figures to disk.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, get_file_names, check_or_create_dir
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=20) 
mpl.rc('axes', labelsize=20) 
mpl.rc('xtick', labelsize=20) 
mpl.rc('ytick', labelsize=20) 
mpl.rc('legend', fontsize=16) 

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
apo_stages = ["early", "mid", "late"]
TIMES = ["48hr", "72hr", "96hr"]

path_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/values/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/figures/"
check_or_create_dir(path_figures)

for EXP in EXPERIMENTS:
    path_results_exp = path_results+"{}/".format(EXP)
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)

    MEAN_DISTS = np.zeros((len(apo_stages), len(TIMES)))
    STDS_DISTS = np.zeros((len(apo_stages), len(TIMES)))

    MEANS = []
    MEANS_APO = []

    MEAN_DISTS_APO = np.zeros((len(apo_stages), len(TIMES)))
    STDS_DISTS_APO = np.zeros((len(apo_stages), len(TIMES)))

    for aps, apo_stage in enumerate(apo_stages):
        means = []
        means_apo = []

        DISTS = []
        DISTS_pre = []

        DISTS_apo = []
        DISTS_apo_pre = []

        DISTS_CONT = []
        DISTS_CONT_apo = []

        for T, TIME in enumerate(TIMES):
            dists_contour_Casp3 = []
            dists_contour_A12 = []
            dists_contour_F3 = []
            dists_centroid_Casp3 = []
            dists_centroid_A12 = []
            dists_centroid_F3 = []

            means.append([])
            means_apo.append([])

            ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
            path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/WT/'.format(EXP, TIME)
            path_save_results='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, "WT")

            ### GET FULL FILE NAME AND FILE CODE ###
            files = get_file_names(path_data_dir)

            if EXP=="2023_11_17_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                if TIME=="96hr":
                    channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
            elif EXP=="2024_03_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]


            for f, file in enumerate(files):
                path_data = path_data_dir+file
                file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)

                file_path = path_save_results+embcode
                dists_contour_Casp3_current = np.load(file_path+"_dists_contour_Casp3.npy")
                dists_contour_A12_current = np.load(file_path+"_dists_contour_A12.npy")
                dists_contour_F3_current = np.load(file_path+"_dists_contour_F3.npy")
                
                dists_centroid_Casp3_current = np.load(file_path+"_dists_centroid_Casp3.npy")
                dists_centroid_A12_current = np.load(file_path+"_dists_centroid_A12.npy")
                dists_centroid_F3_current = np.load(file_path+"_dists_centroid_F3.npy")
                
                dists_contour_Casp3 = [*dists_contour_Casp3, *dists_contour_Casp3_current]
                dists_contour_A12 = [*dists_contour_A12, *dists_contour_A12_current]
                dists_contour_F3 = [*dists_contour_F3, *dists_contour_F3_current]
                
                dists_centroid_Casp3 = [*dists_centroid_Casp3, *dists_centroid_Casp3_current]
                dists_centroid_A12 = [*dists_centroid_A12, *dists_centroid_A12_current]
                dists_centroid_F3 = [*dists_centroid_F3, *dists_centroid_F3_current]

                dists_contour = np.array([*dists_contour_A12, *dists_contour_F3])
                dists_centroid = np.array([*dists_centroid_A12, *dists_centroid_F3])
                dists = dists_centroid / (dists_centroid + dists_contour)
                dists_apo = np.array(dists_centroid_Casp3) / (np.array(dists_centroid_Casp3) + np.array(dists_contour_Casp3))

                means[-1].append(np.mean(dists))
                means_apo[-1].append(np.mean(dists_apo))

            ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
            path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/KO/'.format(EXP, TIME)
            path_save_results='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, "KO")

            ### GET FULL FILE NAME AND FILE CODE ###
            files = get_file_names(path_data_dir)

            if EXP=="2023_11_17_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                if TIME=="96hr":
                    channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
            elif EXP=="2024_03_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

            for f, file in enumerate(files):
                path_data = path_data_dir+file
                file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)

                file_path = path_save_results+embcode
                dists_contour_Casp3_current = np.load(file_path+"_dists_contour_Casp3.npy")
                dists_contour_A12_current = np.load(file_path+"_dists_contour_A12.npy")
                dists_contour_F3_current = np.load(file_path+"_dists_contour_F3.npy")
                
                dists_centroid_Casp3_current = np.load(file_path+"_dists_centroid_Casp3.npy")
                dists_centroid_A12_current = np.load(file_path+"_dists_centroid_A12.npy")
                dists_centroid_F3_current = np.load(file_path+"_dists_centroid_F3.npy")
                
                dists_contour_Casp3 = [*dists_contour_Casp3, *dists_contour_Casp3_current]
                dists_contour_A12 = [*dists_contour_A12, *dists_contour_A12_current]
                dists_contour_F3 = [*dists_contour_F3, *dists_contour_F3_current]
                
                dists_centroid_Casp3 = [*dists_centroid_Casp3, *dists_centroid_Casp3_current]
                dists_centroid_A12 = [*dists_centroid_A12, *dists_centroid_A12_current]
                dists_centroid_F3 = [*dists_centroid_F3, *dists_centroid_F3_current]
                
                dists_contour = np.array([*dists_contour_A12, *dists_contour_F3])
                dists_centroid = np.array([*dists_centroid_A12, *dists_centroid_F3])
                dists = dists_centroid / (dists_centroid + dists_contour)
                dists_apo = np.array(dists_centroid_Casp3) / (np.array(dists_centroid_Casp3) + np.array(dists_contour_Casp3))

                means[-1].append(np.mean(dists))
                means_apo[-1].append(np.mean(dists_apo))

            dists_contour = np.array([*dists_contour_A12, *dists_contour_F3])
            dists_centroid = np.array([*dists_centroid_A12, *dists_centroid_F3])

            dists = dists_centroid / (dists_centroid + dists_contour)
            dists_apo = np.array(dists_centroid_Casp3) / (np.array(dists_centroid_Casp3) + np.array(dists_contour_Casp3))

            DISTS.append(dists)
            DISTS_CONT.append(dists_contour)
            DISTS_apo.append(dists_apo)
            DISTS_CONT_apo.append(dists_contour_Casp3)

        mean_dists = [np.mean(means[t]) for t in range(len(TIMES))]
        mean_dists_apo = [np.mean(means_apo[t]) for t in range(len(TIMES))]

        stds_dists = [np.std(means[t]) for t in range(len(TIMES))]
        stds_dists_apo = [np.std(means_apo[t]) for t in range(len(TIMES))]

        MEAN_DISTS[aps] = np.array(mean_dists) 
        MEAN_DISTS_APO[aps] = np.array(mean_dists_apo) 

        STDS_DISTS[aps] = np.array(stds_dists) 
        STDS_DISTS_APO[aps] = np.array(stds_dists_apo) 

        MEANS.append(means)
        MEANS_APO.append(means_apo)


    fig, ax = plt.subplots(2,3, figsize=(15,8))
    for t, TIME in enumerate(TIMES):
        apo_mean = plt.Circle((0,0), MEAN_DISTS[0, t], linestyle='--', color="k", fill=False, linewidth=3, label=r"$G(x)$")
        ax[0, t].add_patch(apo_mean)
        apo_mean = plt.Circle((0,0), MEAN_DISTS[0, t], linestyle='--', color="k", fill=False, linewidth=3, label=r"$G(x)$")
        ax[1, t].add_patch(apo_mean)
        theta = np.linspace(0, 2 * np.pi, 100)
        
        # Parametric equations for the circles
        r1 = MEAN_DISTS[0, t] - STDS_DISTS[0, t]
        r2 = MEAN_DISTS[0, t] + STDS_DISTS[0, t]

        x1 = 0 + r1 * np.cos(theta)
        y1 = 0 + r1 * np.sin(theta)
        
        x2 = 0 + r2 * np.cos(theta)
        y2 = 0 + r2 * np.sin(theta)
        ax[0, t].fill(np.concatenate([x1, x2[::-1]]), np.concatenate([y1, y2[::-1]]), color="k", alpha=0.2, label=r"$\sigma_M \ G(x)$")
        ax[1, t].fill(np.concatenate([x1, x2[::-1]]), np.concatenate([y1, y2[::-1]]), color="k", alpha=0.2, label=r"$\sigma_M \ G(x)$")

        for aps, apo_stage in enumerate(apo_stages):
            apo_mean = plt.Circle((0,0), MEAN_DISTS_APO[aps, t], color="C{}".format(aps), fill=False, linewidth=3, label=r"$F(x)$ {} apoptosis".format(apo_stage))
            ax[0, t].add_patch(apo_mean)
                # Generate values for the angle theta
            theta = np.linspace(0, 2 * np.pi, 100)
            
            # Parametric equations for the circles
            r1 = MEAN_DISTS_APO[aps, t] - STDS_DISTS_APO[aps, t]
            r2 = MEAN_DISTS_APO[aps, t] + STDS_DISTS_APO[aps, t]

            x1 = 0 + r1 * np.cos(theta)
            y1 = 0 + r1 * np.sin(theta)
            
            x2 = 0 + r2 * np.cos(theta)
            y2 = 0 + r2 * np.sin(theta)
            
            # Fill the area between the circles
            # ax[0, t].fill_between(x1, y2, y1, color="C{}".format(aps), alpha=0.2)
            # Fill the area between the circles
            # ax[0, t].fill(x2, y2, x1, y1, color="C{}".format(aps), alpha=0.5)
            ax[0, t].fill(np.concatenate([x1, x2[::-1]]), np.concatenate([y1, y2[::-1]]), color="C{}".format(aps), alpha=0.2, label=r"$\sigma_M \ F(x)$ {} apoptosis".format(apo_stage))
            apo_mean = plt.Circle((0,0), MEAN_DISTS_APO[aps, t], color="C{}".format(aps), fill=False, linewidth=3, label=apo_stage)
            ax[1, t].add_patch(apo_mean)
            ax[1, t].fill(np.concatenate([x1, x2[::-1]]), np.concatenate([y1, y2[::-1]]), color="C{}".format(aps), alpha=0.2, label=r"$\sigma_M \ G(x)$ {} apoptosis".format(apo_stage))


        gastruloid_border = plt.Circle((0,0), 1.0, color="k", fill=False, linewidth=4, label="Gastruloid edge")
        ax[0, t].add_patch(gastruloid_border)
        ax[0, t].scatter([0], [0], color="k",label="Gatruloid centroid")
        ax[0, t].set_xlim(-1.2, 1.2)
        ax[0, t].set_ylim(-1.2, 1.2)
        ax[0, t].set_title(TIME)
        ax[0, t].set_aspect('equal')
        ax[0, t].grid(False)
        ax[0, t].axis('off')

        # inset

        xlims = np.array([0.2, 0.9])
        ylims = np.array([-0.9, -0.2])
        rect = plt.Rectangle((xlims[0], ylims[0]), np.diff(xlims), np.diff(ylims), color="grey", fill=False, linewidth=4)
        ax[0, t].add_patch(rect)

        rect = plt.Rectangle((xlims[0], ylims[0]), np.diff(xlims), np.diff(ylims), color="grey", fill=False, linewidth=4)
        ax[1, t].add_patch(rect)

        xy1 = (xlims[0],ylims[1])
        xy2 = (xlims[0],ylims[0])

        con = ConnectionPatch(xyA=xy1, xyB=xy2, coordsA="data", coordsB="data",
                            axesA=ax[1, t], axesB=ax[0, t], color="grey", linewidth=4)
        ax[1, t].add_artist(con)

        xy2 = (xlims[1],ylims[0])
        xy1 = (xlims[1],ylims[1])

        con = ConnectionPatch(xyA=xy1, xyB=xy2, coordsA="data", coordsB="data",
                            axesA=ax[1, t], axesB=ax[0, t], color="grey", linewidth=4)
        ax[1, t].add_artist(con)

        gastruloid_border = plt.Circle((0,0), 1.0, color="k", fill=False, linewidth=4, label="gastruloid edge")
        ax[1, t].add_patch(gastruloid_border)

        ax[1, t].scatter([0], [0], color="k",label="gatruloid centroid")
        ax[1, t].set_xlim((xlims[0] - 0.005, xlims[1] + 0.005))
        ax[1, t].set_ylim((ylims[0] - 0.005, ylims[1] + 0.005))
        # ax[1, t].set_title(TIME)
        ax[1, t].set_aspect('equal')
        ax[1, t].grid(False)
        ax[1, t].axis('off')

    ax[0,-1].legend(bbox_to_anchor=(1.1, 1.0))
    plt.tight_layout()
    plt.savefig(path_figures_exp+"radialevo.svg")
    plt.savefig(path_figures_exp+"radialevo.pdf")
    plt.show()