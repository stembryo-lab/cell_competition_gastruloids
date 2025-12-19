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
apo_stage = apo_stages[0]
TIMES = ["48hr", "72hr", "96hr"]

path_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/values/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/figures/"
check_or_create_dir(path_figures)

DISTS_F3_WT = []
DISTS_A12_WT = []
DISTS_apo_WT = []

DISTS_F3_KO = []
DISTS_A12_KO = []
DISTS_apo_KO = []

all_files = []
all_data = []

for EXP in EXPERIMENTS:
    path_results_exp = path_results+"{}/".format(EXP)
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)

    for TIME in TIMES:
        dists_F3 = []
        dists_A12 = []
        dists_Casp3 = []
        
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
            
            dists = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            all_files.append(file + " F3")
            all_data.append(dists)
            dists_F3 = [*dists_F3, *dists]
            
            dists = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            all_files.append(file + " A12")
            all_data.append(dists)
            dists_A12 = [*dists_A12, *dists]
            
            dists = dists_centroid_Casp3_current / (dists_centroid_Casp3_current + dists_contour_Casp3_current)
            dists_Casp3 = [*dists_Casp3, *dists]
            
                    
        DISTS_F3_KO.append(dists_F3)
        DISTS_A12_KO.append(dists_A12)
        DISTS_apo_KO.append(dists_Casp3)

        ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
        path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/KO/'.format(EXP, TIME)
        path_save_results='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, "KO")
        
        ### GET FULL FILE NAME AND FILE CODE ###
        files = get_file_names(path_data_dir)

        dists_F3 = []
        dists_A12 = []
        dists_Casp3 = []
        
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
            
            dists = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            all_files.append(file + " F3")
            all_data.append(dists)
            dists_F3 = [*dists_F3, *dists]
            
            dists = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            all_files.append(file + " A12")
            all_data.append(dists)
            dists_A12 = [*dists_A12, *dists]
            
            dists = dists_centroid_Casp3_current / (dists_centroid_Casp3_current + dists_contour_Casp3_current)
            dists_Casp3 = [*dists_Casp3, *dists]
        

        DISTS_F3_WT.append(dists_F3)
        DISTS_A12_WT.append(dists_A12)
        DISTS_apo_WT.append(dists_Casp3)

    density=True
    fig, ax = plt.subplots(2,3, figsize=(12,6),sharex=True)
    for t, TIME in enumerate(TIMES):
        
        ax[0,t].hist(DISTS_F3_WT[t], color="green", alpha=0.5, bins=50, density=density)
        ax[0,t].hist(DISTS_A12_WT[t], color="magenta", alpha=0.5, bins=50, density=density)
        ax[0,t].set_yticks([])
        # ax[0,t].set_xlim(-0.1,1.1)
        ax[0,t].set_title(TIME)
        ax[0,t].spines[['left', 'right', 'top']].set_visible(False)
        if t==0:
            ax[0,t].set_ylabel("WT")

        ax[1,t].hist(DISTS_F3_KO[t], color="green", alpha=0.5, bins=50, density=density, label="F3")
        ax[1,t].hist(DISTS_A12_KO[t], color="magenta", alpha=0.5, bins=50, density=density, label="A12")
        ax[1,t].set_yticks([])
        # ax[1,t].set_xlim(-0.1,1.1)
        ax[1,t].spines[['left', 'right', 'top']].set_visible(False)
        ax[1,t].set_xlabel(r"relative position on gastruloid")

        if t ==0:
            ax[1,t].set_ylabel("KO")

        if t==len(TIMES)-1:
            ax[1,t].legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(path_figures_exp+"dists_condition.svg")
    plt.savefig(path_figures_exp+"dists_condition.pdf")
    plt.show()
        
    density=True
    bins=20
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,2, figsize=(10,3.5))
    for t, TIME in enumerate(TIMES):

        dists = [*DISTS_F3_WT[t], *DISTS_A12_WT[t], *DISTS_F3_KO[t], *DISTS_A12_KO[t]]
        ax[0].hist(dists, alpha=0.5, bins=bins, density=density, label=TIME)
        ax[0].set_yticks([])
        ax[0].set_xlim(-0.1,1.1)
        ax[0].set_xlabel(r"relative position on gastruloid ($P$)")
        ax[0].spines[['left', 'right', 'top']].set_visible(False)
        ax[0].set_title(r"$G(x)$")

        dists = [*DISTS_apo_WT[t], *DISTS_apo_KO[t]]
        ax[1].hist(dists, alpha=0.5, bins=bins, density=density, label=TIME)
        ax[1].set_yticks([])
        ax[1].set_xlim(-0.1,1.1)
        ax[1].set_xlabel(r"relative position on gastruloid ($P$)")
        ax[1].spines[['left', 'right', 'top']].set_visible(False)
        ax[1].set_title(r"$F(x)$")
        if t==len(TIMES)-1:
            ax[1].legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(path_figures_exp+"dists_apo.svg")
    plt.savefig(path_figures_exp+"dists_apo.pdf")
    plt.show()
