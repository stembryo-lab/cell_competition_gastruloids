"""
Binned radial distribution analysis.

Loads per-cell radial-position values (P = d_centroid / (d_centroid + d_edge)),
bins them into radial shells (multiple bin numbers), and exports:
- per-embryo binned profiles to CSV
- summary plots (WT vs KO; F3 vs A12) in both raw and volume-normalized form.
Also writes single-cell radial-position tables for each condition/time point.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, get_file_names, check_or_create_dir
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

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
path_csvs = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/csvs/"
check_or_create_dir(path_results)
check_or_create_dir(path_figures)
check_or_create_dir(path_csvs)


for EXP in EXPERIMENTS:
    path_results_exp = path_results+"{}/".format(EXP)
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)
    path_csvs_exp = path_csvs+"{}/".format(EXP)
    check_or_create_dir(path_csvs_exp)

    DISTS_F3_WT = []
    DISTS_A12_WT = []
    DISTS_F3_apo_WT = []
    DISTS_A12_apo_WT = []

    DISTS_F3_KO = []
    DISTS_A12_KO = []
    DISTS_F3_apo_KO = []
    DISTS_A12_apo_KO = []

    path_figures_exp_norm = path_figures_exp+"norm/"
    check_or_create_dir(path_figures_exp_norm)

    path_figures_exp_raw = path_figures_exp+"raw/"
    check_or_create_dir(path_figures_exp_raw)

    path_csvs_exp_norm = path_csvs_exp+"norm/"
    check_or_create_dir(path_csvs_exp_norm)

    path_csvs_exp_raw = path_csvs_exp+"raw/"
    check_or_create_dir(path_csvs_exp_raw)
    
    if EXP=="2023_11_17_Casp3":
        TIMES = ["48hr", "72hr", "96hr"]
    else:        
        TIMES = ["48hr", "72hr", "96hr"]

    g_files_WT = []
    g_files_KO = []

    for t, TIME in enumerate(TIMES):
        dists_F3 = []
        dists_A12 = []
        
        ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
        path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/KO/'.format(EXP, TIME)
        path_save_results_early='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, "KO")
        files = get_file_names(path_data_dir)

        g_files_KO.append([])
        g_files_WT.append([])
        
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if "96hr" in path_data_dir:
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        else:
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

        for f, file in enumerate(files):
            path_data = path_data_dir+file
            file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
            
            file_path_early = path_save_results_early+embcode
            
            dists_contour_A12_current = np.load(file_path_early+"_dists_contour_A12.npy")
            dists_contour_F3_current = np.load(file_path_early+"_dists_contour_F3.npy")
            
            dists_centroid_A12_current = np.load(file_path_early+"_dists_centroid_A12.npy")
            dists_centroid_F3_current = np.load(file_path_early+"_dists_centroid_F3.npy")
            
            dists_f3 = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            dists_F3.append(dists_f3)
            
            dists_a12 = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            dists_A12.append(dists_a12)
            
            g_files_KO[t].append(embcode)
            
        DISTS_F3_KO.append(dists_F3)
        DISTS_A12_KO.append(dists_A12)
        
        ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
        path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/WT/'.format(EXP, TIME)
        path_save_results_early='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, "WT")
        files = get_file_names(path_data_dir)

        dists_F3 = []
        dists_A12 = []
        dists_Casp3 = []
        
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if "96hr" in path_data_dir:
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        else:
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

        for f, file in enumerate(files):
            path_data = path_data_dir+file
            file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)

            file_path_early = path_save_results_early+embcode
            
            dists_contour_A12_current = np.load(file_path_early+"_dists_contour_A12.npy")
            dists_contour_F3_current = np.load(file_path_early+"_dists_contour_F3.npy")

            dists_centroid_A12_current = np.load(file_path_early+"_dists_centroid_A12.npy")
            dists_centroid_F3_current = np.load(file_path_early+"_dists_centroid_F3.npy")
            
            dists_f3 = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            dists_F3.append(dists_f3)
            
            dists_a12 = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            dists_A12.append(dists_a12)
            
            g_files_WT[t].append(embcode)
            
        DISTS_F3_WT.append(dists_F3)
        DISTS_A12_WT.append(dists_A12)

    sct_size=10
    line_w = 2

    t = 0
    bin_step = 0.1
    bins = np.arange(0, 1.0+bin_step, bin_step)

    total_bins_all = [6, 11, 16, 21]
    for total_bins in total_bins_all:
        path_figures_bins = path_figures_exp_norm+"{}/".format(total_bins)
        check_or_create_dir(path_figures_bins)
   
        path_csvs_bins = path_csvs_exp_norm+"{}/".format(total_bins)
        check_or_create_dir(path_csvs_bins)
             
        bins= np.linspace(0, 1.0, num=total_bins)
        x = bins[:-1] + (np.diff(bins)/2)
        fig, ax = plt.subplots(2,3, figsize=(14,8),sharex=True, sharey=False)
        for t, TIME in enumerate(TIMES):
            path_figures_time = path_figures_bins+"{}/".format(TIME)
            check_or_create_dir(path_figures_time)

            path_csvs_time = path_csvs_bins+"{}/".format(TIME)
            check_or_create_dir(path_csvs_time)
                    
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_F3_WT[t])):
                _counts, bins = np.histogram(DISTS_F3_WT[t][g], bins=bins, density=True)
                counts = [_counts[i]/((4/3)*np.pi*(bins[i+1]**3-bins[i]**3)) for i in range(len(_counts))]
                ax[0,t].scatter(x, counts, color="green", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            f_names = g_files_WT[t]
            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"F3_WT.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[0,t].plot(x, counts_mean, color="green", lw=line_w)
            ax[0,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="green", alpha=0.5)

            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_A12_WT[t])):
                _counts, bins = np.histogram(DISTS_A12_WT[t][g], bins=bins, density=True)
                counts = [_counts[i]/((4/3)*np.pi*(bins[i+1]**3-bins[i]**3)) for i in range(len(_counts))]
                ax[0,t].scatter(x, counts, color="magenta", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"A12_WT.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[0,t].plot(x, counts_mean, color="magenta", lw=line_w)
            ax[0,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="magenta", alpha=0.5)
            
            ax[0,t].set_xlim(-0.1,1.1)
            ax[0,t].set_title(TIME)
            
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_F3_KO[t])):
                _counts, bins = np.histogram(DISTS_F3_KO[t][g], bins=bins, density=True)
                counts = [_counts[i]/((4/3)*np.pi*(bins[i+1]**3-bins[i]**3)) for i in range(len(_counts))]
                ax[1,t].scatter(x, counts, color="green", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            f_names = g_files_KO[t]
            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"F3_KO.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[1,t].plot(x, counts_mean, color="green", lw=line_w)
            ax[1,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="green", alpha=0.5)
            
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_A12_KO[t])):
                _counts, bins = np.histogram(DISTS_A12_KO[t][g], bins=bins, density=True)
                counts = [_counts[i]/((4/3)*np.pi*(bins[i+1]**3-bins[i]**3)) for i in range(len(_counts))]
                ax[1,t].scatter(x, counts, color="magenta", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"A12_KO.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[1,t].plot(x, counts_mean, color="magenta", lw=line_w)
            ax[1,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="magenta", alpha=0.5)
            
            ax[1,t].set_xlim(-0.1,1.1)
            # ax[1,t].spines[['left', 'right', 'top']].set_visible(False)
            ax[1,t].set_xlabel(r"relative position on gastruloid")
            if t==0:
                ax[0,t].set_ylabel("WT")
                ax[1,t].set_ylabel("KO")

        plt.tight_layout()
        plt.savefig(path_figures_bins+"radial_dist_binned.svg")

    total_bins_all = [6, 11, 16, 21]
    for total_bins in total_bins_all:
        path_figures_bins = path_figures_exp_raw+"{}/".format(total_bins)
        check_or_create_dir(path_figures_bins)

        path_csvs_bins = path_csvs_exp_raw+"{}/".format(total_bins)
        check_or_create_dir(path_csvs_bins)
         
        bins= np.linspace(0, 1.0, num=total_bins)
        x = bins[:-1] + (np.diff(bins)/2)
        fig, ax = plt.subplots(2,3, figsize=(14,8),sharex=True, sharey=False)
        for t, TIME in enumerate(TIMES):
            path_figures_time = path_figures_bins+"{}/".format(TIME)
            check_or_create_dir(path_figures_time)
            
            path_csvs_time = path_csvs_bins+"{}/".format(TIME)
            check_or_create_dir(path_csvs_time)
        
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_F3_WT[t])):
                counts, bins = np.histogram(DISTS_F3_WT[t][g], bins=bins, density=True)
                ax[0,t].scatter(x, counts, color="green", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            f_names = g_files_WT[t]
            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"F3_WT.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[0,t].plot(x, counts_mean, color="green", lw=line_w)
            ax[0,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="green", alpha=0.5)

            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_A12_WT[t])):
                counts, bins = np.histogram(DISTS_A12_WT[t][g], bins=bins, density=True)
                ax[0,t].scatter(x, counts, color="magenta", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"A12_WT.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[0,t].plot(x, counts_mean, color="magenta", lw=line_w)
            ax[0,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="magenta", alpha=0.5)
            
            ax[0,t].set_xlim(-0.1,1.1)
            ax[0,t].set_title(TIME)
            
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_F3_KO[t])):
                counts, bins = np.histogram(DISTS_F3_KO[t][g], bins=bins, density=True)
                ax[1,t].scatter(x, counts, color="green", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            f_names = g_files_KO[t]
            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"F3_KO.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[1,t].plot(x, counts_mean, color="green", lw=line_w)
            ax[1,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="green", alpha=0.5)
            
            COUNTS = []
            x = bins[:-1] + (np.diff(bins)/2)
            for g in range(len(DISTS_A12_KO[t])):
                counts, bins = np.histogram(DISTS_A12_KO[t][g], bins=bins, density=True)
                ax[1,t].scatter(x, counts, color="magenta", s=sct_size, alpha=0.5)
                COUNTS.append(counts)

            df = pd.DataFrame(COUNTS, columns=x, index=f_names)
            df.to_csv(path_csvs_time+"A12_KO.csv")
            
            counts_mean = np.mean(COUNTS, axis=0)
            counts_std = np.std(COUNTS, axis=0)
            ax[1,t].plot(x, counts_mean, color="magenta", lw=line_w)
            ax[1,t].fill_between(x, counts_mean-counts_std, counts_mean+counts_std, color="magenta", alpha=0.5)
            
            ax[1,t].set_xlim(-0.1,1.1)
            # ax[1,t].spines[['left', 'right', 'top']].set_visible(False)
            ax[1,t].set_xlabel(r"relative position on gastruloid")
            if t==0:
                ax[0,t].set_ylabel("WT")
                ax[1,t].set_ylabel("KO")

        plt.tight_layout()
        plt.savefig(path_figures_bins+"radial_dist_binned.svg")

    path_csvs_single_cell = path_csvs+"single_cell/"
    check_or_create_dir(path_csvs_single_cell)
     
    for TIME in TIMES:
        path_csvs_time = path_csvs_single_cell+"{}/".format(TIME)
        check_or_create_dir(path_csvs_time)
        
        path_csvs_KO = path_csvs_time+"{}/".format("KO")
        check_or_create_dir(path_csvs_KO)
        
        path_csvs_WT = path_csvs_time+"{}/".format("WT")
        check_or_create_dir(path_csvs_WT)
        
        ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
        path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/KO/'.format(EXP, TIME)
        path_save_results_early='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[0], TIME, "KO")
        path_save_results_mid='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[1], TIME, "KO")
        path_save_results_late='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[2], TIME, "KO")
        files = get_file_names(path_data_dir)
        
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if "96hr" in path_data_dir:
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        else:
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

        for f, file in enumerate(files):
            path_data = path_data_dir+file
            file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
            
            file_path_early = path_save_results_early+embcode
            dists_contour_Casp3_A12_current_early = np.load(file_path_early+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_early = np.load(file_path_early+"_dists_contour_Casp3_F3.npy")
            
            file_path_mid = path_save_results_mid+embcode
            dists_contour_Casp3_A12_current_mid = np.load(file_path_mid+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_mid = np.load(file_path_mid+"_dists_contour_Casp3_F3.npy")
            
            file_path_late = path_save_results_late+embcode
            dists_contour_Casp3_A12_current_late = np.load(file_path_late+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_late = np.load(file_path_late+"_dists_contour_Casp3_F3.npy")


            dists_centroid_Casp3_A12_current_early = np.load(file_path_early+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_early = np.load(file_path_early+"_dists_centroid_Casp3_F3.npy")
            
            file_path_mid = path_save_results_mid+embcode
            dists_centroid_Casp3_A12_current_mid = np.load(file_path_mid+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_mid = np.load(file_path_mid+"_dists_centroid_Casp3_F3.npy")
            
            file_path_late = path_save_results_late+embcode
            dists_centroid_Casp3_A12_current_late = np.load(file_path_late+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_late = np.load(file_path_late+"_dists_centroid_Casp3_F3.npy")
            
            dists_contour_A12_current = np.load(file_path_early+"_dists_contour_A12.npy")
            dists_contour_F3_current = np.load(file_path_early+"_dists_contour_F3.npy")

            dists_centroid_A12_current = np.load(file_path_early+"_dists_centroid_A12.npy")
            dists_centroid_F3_current = np.load(file_path_early+"_dists_centroid_F3.npy")
            
            dists_f3 = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            dists_a12 = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            
            dists_casp3_f3_early = dists_centroid_Casp3_F3_current_early / (dists_centroid_Casp3_F3_current_early + dists_contour_Casp3_F3_current_early)
            dists_casp3_a12_early = dists_centroid_Casp3_A12_current_early / (dists_centroid_Casp3_A12_current_early + dists_contour_Casp3_A12_current_early)
            
            dists_casp3_f3_mid = dists_centroid_Casp3_F3_current_mid / (dists_centroid_Casp3_F3_current_mid + dists_contour_Casp3_F3_current_mid)
            dists_casp3_a12_mid = dists_centroid_Casp3_A12_current_mid / (dists_centroid_Casp3_A12_current_mid + dists_contour_Casp3_A12_current_mid)
            
            dists_casp3_f3_late = dists_centroid_Casp3_F3_current_late / (dists_centroid_Casp3_F3_current_late + dists_contour_Casp3_F3_current_late)
            dists_casp3_a12_late = dists_centroid_Casp3_A12_current_late / (dists_centroid_Casp3_A12_current_late + dists_contour_Casp3_A12_current_late)
            
            max_len = max(len(dists_f3), len(dists_a12), len(dists_casp3_f3_early), len(dists_casp3_a12_early), len(dists_casp3_f3_mid), len(dists_casp3_a12_mid), len(dists_casp3_f3_late), len(dists_casp3_a12_late))

            F3_pad = list(dists_f3) + [""] * (max_len - len(dists_f3))
            A12_pad = list(dists_a12) + [""] * (max_len - len(dists_a12))
            apo_F3_pad_early = list(dists_casp3_f3_early) + [""] * (max_len - len(dists_casp3_f3_early))
            apo_A12_pad_early = list(dists_casp3_a12_early) + [""] * (max_len - len(dists_casp3_a12_early))

            apo_F3_pad_mid = list(dists_casp3_f3_mid) + [""] * (max_len - len(dists_casp3_f3_mid))
            apo_A12_pad_mid = list(dists_casp3_a12_mid) + [""] * (max_len - len(dists_casp3_a12_mid))

            apo_F3_pad_late = list(dists_casp3_f3_late) + [""] * (max_len - len(dists_casp3_f3_late))
            apo_A12_pad_late = list(dists_casp3_a12_late) + [""] * (max_len - len(dists_casp3_a12_late))

            df = pd.DataFrame({
                "F3": F3_pad,
                "A12": A12_pad,
                "early_apo_F3": apo_F3_pad_early,
                "early_apo_A12": apo_A12_pad_early,
                "mid_apo_F3": apo_F3_pad_mid,
                "mid_apo_A12": apo_A12_pad_mid,
                "late_apo_F3": apo_F3_pad_late,
                "late_apo_A12": apo_A12_pad_late
            })

            df.to_csv(path_csvs_KO+"raw_data.csv", index=False)
        
        ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
        path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/WT/'.format(EXP, TIME)
        path_save_results_early='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[0], TIME, "WT")
        path_save_results_mid='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[1], TIME, "WT")
        path_save_results_late='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stages[2], TIME, "WT")
        files = get_file_names(path_data_dir)
        
        ### GET FULL FILE NAME AND FILE CODE ###
        files = get_file_names(path_data_dir)

        dists_F3 = []
        dists_A12 = []
        dists_Casp3 = []

        for f, file in enumerate(files):
            path_data = path_data_dir+file
            file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)

            file_path_early = path_save_results_early+embcode
            dists_contour_Casp3_A12_current_early = np.load(file_path_early+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_early = np.load(file_path_early+"_dists_contour_Casp3_F3.npy")
            
            file_path_mid = path_save_results_mid+embcode
            dists_contour_Casp3_A12_current_mid = np.load(file_path_mid+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_mid = np.load(file_path_mid+"_dists_contour_Casp3_F3.npy")
            
            file_path_late = path_save_results_late+embcode
            dists_contour_Casp3_A12_current_late = np.load(file_path_late+"_dists_contour_Casp3_A12.npy")
            dists_contour_Casp3_F3_current_late = np.load(file_path_late+"_dists_contour_Casp3_F3.npy")


            dists_centroid_Casp3_A12_current_early = np.load(file_path_early+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_early = np.load(file_path_early+"_dists_centroid_Casp3_F3.npy")
            
            file_path_mid = path_save_results_mid+embcode
            dists_centroid_Casp3_A12_current_mid = np.load(file_path_mid+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_mid = np.load(file_path_mid+"_dists_centroid_Casp3_F3.npy")
            
            file_path_late = path_save_results_late+embcode
            dists_centroid_Casp3_A12_current_late = np.load(file_path_late+"_dists_centroid_Casp3_A12.npy")
            dists_centroid_Casp3_F3_current_late = np.load(file_path_late+"_dists_centroid_Casp3_F3.npy")
            
            dists_contour_A12_current = np.load(file_path_early+"_dists_contour_A12.npy")
            dists_contour_F3_current = np.load(file_path_early+"_dists_contour_F3.npy")

            dists_centroid_A12_current = np.load(file_path_early+"_dists_centroid_A12.npy")
            dists_centroid_F3_current = np.load(file_path_early+"_dists_centroid_F3.npy")
            
            dists_f3 = dists_centroid_F3_current / (dists_centroid_F3_current + dists_contour_F3_current)
            dists_a12 = dists_centroid_A12_current / (dists_centroid_A12_current + dists_contour_A12_current)
            
            dists_casp3_f3_early = dists_centroid_Casp3_F3_current_early / (dists_centroid_Casp3_F3_current_early + dists_contour_Casp3_F3_current_early)
            dists_casp3_a12_early = dists_centroid_Casp3_A12_current_early / (dists_centroid_Casp3_A12_current_early + dists_contour_Casp3_A12_current_early)
            
            dists_casp3_f3_mid = dists_centroid_Casp3_F3_current_mid / (dists_centroid_Casp3_F3_current_mid + dists_contour_Casp3_F3_current_mid)
            dists_casp3_a12_mid = dists_centroid_Casp3_A12_current_mid / (dists_centroid_Casp3_A12_current_mid + dists_contour_Casp3_A12_current_mid)
            
            dists_casp3_f3_late = dists_centroid_Casp3_F3_current_late / (dists_centroid_Casp3_F3_current_late + dists_contour_Casp3_F3_current_late)
            dists_casp3_a12_late = dists_centroid_Casp3_A12_current_late / (dists_centroid_Casp3_A12_current_late + dists_contour_Casp3_A12_current_late)
            
            max_len = max(len(dists_f3), len(dists_a12), len(dists_casp3_f3_early), len(dists_casp3_a12_early), len(dists_casp3_f3_mid), len(dists_casp3_a12_mid), len(dists_casp3_f3_late), len(dists_casp3_a12_late))

            F3_pad = list(dists_f3) + [""] * (max_len - len(dists_f3))
            A12_pad = list(dists_a12) + [""] * (max_len - len(dists_a12))
            apo_F3_pad_early = list(dists_casp3_f3_early) + [""] * (max_len - len(dists_casp3_f3_early))
            apo_A12_pad_early = list(dists_casp3_a12_early) + [""] * (max_len - len(dists_casp3_a12_early))

            apo_F3_pad_mid = list(dists_casp3_f3_mid) + [""] * (max_len - len(dists_casp3_f3_mid))
            apo_A12_pad_mid = list(dists_casp3_a12_mid) + [""] * (max_len - len(dists_casp3_a12_mid))

            apo_F3_pad_late = list(dists_casp3_f3_late) + [""] * (max_len - len(dists_casp3_f3_late))
            apo_A12_pad_late = list(dists_casp3_a12_late) + [""] * (max_len - len(dists_casp3_a12_late))

            df = pd.DataFrame({
                "F3": F3_pad,
                "A12": A12_pad,
                "early_apo_F3": apo_F3_pad_early,
                "early_apo_A12": apo_A12_pad_early,
                "mid_apo_F3": apo_F3_pad_mid,
                "mid_apo_A12": apo_A12_pad_mid,
                "late_apo_F3": apo_F3_pad_late,
                "late_apo_A12": apo_A12_pad_late
            })

            df.to_csv(path_csvs_WT+"raw_data.csv", index=False)
