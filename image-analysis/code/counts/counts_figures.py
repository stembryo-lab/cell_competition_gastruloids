"""
Compute cell counts and apoptotic proportions from saved segmentations.

Loads the precomputed F3/A12 segmentations and Casp3 segmentations (early/mid/late),
counts objects per embryo, assigns apoptotic events to F3 vs A12 by channel intensity,
and exports resulting plots showing the results.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.rcParams.update({
    "text.usetex": True,
})
mpl.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath'
mpl.rc('font', size=16) 
mpl.rc('axes', labelsize=16) 
mpl.rc('xtick', labelsize=16) 
mpl.rc('ytick', labelsize=16) 
mpl.rc('legend', fontsize=16) 
    

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
CONDS = ["WT", "KO"]
TIMES = ["48hr", "72hr", "96hr"]

APO_STAGES = ["early","mid", "late"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
master_path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/counts/"
check_or_create_dir(master_path_to_save_results)

for EXP in EXPERIMENTS:    
    path_figures_exp = master_path_to_save_results+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)
    
    F3_WT = []
    A12_WT = []
    F3_KO = []
    A12_KO = []
    
    Casp3_F3_KO = {}
    Casp3_A12_KO = {}
    Casp3_F3_WT = {}
    Casp3_A12_WT = {}
    
    for apo_stage in APO_STAGES:
        Casp3_F3_KO[apo_stage] = []
        Casp3_A12_KO[apo_stage] = []
        Casp3_F3_WT[apo_stage] = []
        Casp3_A12_WT[apo_stage] = []
        
    for TIME in TIMES:
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if TIME=="96hr":
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        elif EXP=="2024_03_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
    
        for COND in CONDS:
            path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
            path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
            files = get_file_names(path_data_dir)

            F3_counts = []
            A12_counts = []
            
            Casp3_A12_counts = {}
            Casp3_F3_counts = {}
            for apo_stage in APO_STAGES:
                Casp3_A12_counts[apo_stage] = []
                Casp3_F3_counts[apo_stage] = []
            

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
                
                # total F3 and A12s
                F3_counts.append(len(CT_F3.jitcells))
                A12_counts.append(len(CT_A12.jitcells))
                
                ### CORRECT DISTRIBUTIONS ###

                F3_dist = []
                for cell in CT_F3.jitcells: 
                    for zid, z in enumerate(cell.zs[0]):
                        mask = cell.masks[0][zid]
                        img = CT_F3.hyperstack[0,z, channel_names.index("F3")]
                        F3_dist.append(np.mean(img[mask[:,1], mask[:,0]]))
                        
                A12_dist = []
                for cell in CT_A12.jitcells: 
                    for zid, z in enumerate(cell.zs[0]):
                        mask = cell.masks[0][zid]
                        img = CT_A12.hyperstack[0,z, channel_names.index("A12")]
                        A12_dist.append(np.mean(img[mask[:,1], mask[:,0]]))

                ## Now contract the shape as much as we want. 
                F3_dist = np.array(F3_dist)
                A12_dist = np.array(A12_dist)

                mdiff = np.mean(F3_dist) - np.mean(A12_dist)
                if mdiff > 0:
                    A12_dist += mdiff
                else: 
                    F3_dist -= mdiff 

                for apo_stage in APO_STAGES:

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
                    
                    # Now we have corrected for diferences on the mCherry and the emiRFP channels
                    zres = CT_F3.metadata["Zresolution"]
                    xyres = CT_F3.metadata["XYresolution"]

                    Casp3_F3 = 0
                    Casp3_A12 = 0
                    for cell in CT_Casp3.jitcells:
                        casp3 = []
                        a12 = []
                        f3 = []
                        for zid, z in enumerate(cell.zs[0]):
                            mask = cell.masks[0][zid]
                            img = CT_Casp3.hyperstack[0,z, channel_names.index("Casp3")]
                            casp3.append(np.mean(img[mask[:,1], mask[:,0]]))
                            
                            img = CT_Casp3.hyperstack[0,z, channel_names.index("A12")]
                            a12.append(np.mean(img[mask[:,1], mask[:,0]]))
                            
                            img = CT_Casp3.hyperstack[0,z, channel_names.index("F3")]
                            f3.append(np.mean(img[mask[:,1], mask[:,0]]))

                        zz = np.int64(cell.centers[0][0])
                        idx = cell.zs[0].index(zz)
                        if f3[idx] > a12[idx]:
                            Casp3_F3+=1
                        else:
                            Casp3_A12+=1
                    
                    Casp3_F3_counts[apo_stage].append(Casp3_F3)
                    Casp3_A12_counts[apo_stage].append(Casp3_A12)

            if COND == "WT":
                F3_WT.append(F3_counts)
                A12_WT.append(A12_counts)
                for apo_stage in APO_STAGES:
                    Casp3_F3_WT[apo_stage].append(Casp3_F3_counts[apo_stage])
                    Casp3_A12_WT[apo_stage].append(Casp3_A12_counts[apo_stage])

            elif COND == "KO":
                F3_KO.append(F3_counts)
                A12_KO.append(A12_counts)
                for apo_stage in APO_STAGES:
                    Casp3_F3_KO[apo_stage].append(Casp3_F3_counts[apo_stage])
                    Casp3_A12_KO[apo_stage].append(Casp3_A12_counts[apo_stage])


    path_figure = "{}/counts.pdf".format(path_figures_exp)

    fig, ax = plt.subplots(2,4, figsize=(3.5*4,6))

    ax[0,0].set_title("total cells")
    ax[0,1].set_title("early apoptotis")
    ax[0,2].set_title("mid apoptotis")
    ax[0,3].set_title("late apoptotis")


    ax[0,0].set_ylabel(r"WT \\ cell number")
    ax[1,0].set_ylabel(r"KO \\ cell number")

    ax[0,1].set_ylabel("apoptosis proportion")
    ax[1,1].set_ylabel("apoptosis proportion")

    ### WT ###

    time=range(len(TIMES))

    F3_WT_means = np.array([np.mean(data) for data in F3_WT])
    F3_WT_stds = np.array([np.std(data) for data in F3_WT])

    ax[0,0].plot(time, F3_WT_means, color=[0.0,0.8,0.0], lw=2, zorder=1)
    ax[0,0].scatter(time, F3_WT_means, color=[0.0,0.6,0.0], edgecolor='k', s=50, zorder=2)
    ax[0,0].fill_between(time, F3_WT_means-F3_WT_stds, F3_WT_means+F3_WT_stds, alpha=0.1, color=[0.0,0.8,0.0])

    A12_WT_means = np.array([np.mean(data) for data in A12_WT])
    A12_WT_stds = np.array([np.std(data) for data in A12_WT])

    ax[0,0].plot(time, A12_WT_means, color=[0.8,0.0,0.8], lw=2, zorder=1)
    ax[0,0].scatter(time, A12_WT_means, color=[0.6,0.0,0.6], edgecolor='k', s=50, zorder=2)
    ax[0,0].fill_between(time, A12_WT_means-A12_WT_stds, A12_WT_means+A12_WT_stds, alpha=0.1, color=[0.8,0.0,0.8])


    F3_KO_means = np.array([np.mean(data) for data in F3_KO])
    F3_KO_stds = np.array([np.std(data) for data in F3_KO])

    ax[1,0].plot(time, F3_KO_means, color=[0.0,0.8,0.0], lw=2, zorder=1)
    ax[1,0].scatter(time, F3_KO_means, color=[0.0,0.6,0.0], edgecolor='k', s=50, zorder=2)
    ax[1,0].fill_between(time, F3_KO_means-F3_KO_stds, F3_KO_means+F3_KO_stds, alpha=0.1, color=[0.8,0.0,0.8])

    A12_KO_means = np.array([np.mean(data) for data in A12_KO])
    A12_KO_stds = np.array([np.std(data) for data in A12_KO])

    ax[1,0].plot(time, A12_KO_means, color=[0.8,0.0,0.8], lw=2, zorder=1)
    ax[1,0].scatter(time, A12_KO_means, color=[0.6,0.0,0.6], edgecolor='k', s=50, zorder=2)
    ax[1,0].fill_between(time, A12_KO_means-A12_KO_stds, A12_KO_means+A12_KO_stds, alpha=0.1, color=[0.8,0.0,0.8])
    ax[1,0].sharey(ax[0,0])

    ax[0, 0].set_xticks([i for i in time])
    ax[0, 0].set_xticklabels([None]*len(time))
    ax[1, 0].set_xticks([i for i in time])
    ax[1, 0].set_xticklabels(TIMES)

    ax[0, 0].spines[['right', 'top']].set_visible(False)
    ax[1, 0].spines[['right', 'top']].set_visible(False)

    for st, stage in enumerate(APO_STAGES):

        DATA = [np.array(Casp3_F3_WT[stage][n])/(np.array(F3_WT[n]) + np.array(Casp3_F3_WT[stage][n])) for n in range(len(F3_WT))]
        Casp3_F3_WT_means = np.array([np.mean(data) for data in DATA])
        Casp3_F3_WT_stds = np.array([np.std(data) for data in DATA])

        ax[0,st+1].plot(time, Casp3_F3_WT_means, color=[0.0,0.8,0.0], lw=2, zorder=1)
        ax[0,st+1].scatter(time, Casp3_F3_WT_means, color=[0.0,0.6,0.0], edgecolor='k', s=50, zorder=2)
        ax[0,st+1].fill_between(time, Casp3_F3_WT_means-Casp3_F3_WT_stds, Casp3_F3_WT_means+Casp3_F3_WT_stds, alpha=0.1, color=[0.0,0.8,0.0])

        DATA = [np.array(Casp3_F3_KO[stage][n])/(np.array(F3_KO[n]) + np.array(Casp3_F3_KO[stage][n])) for n in range(len(F3_KO))]
        Casp3_F3_KO_means = np.array([np.mean(data) for data in DATA])
        Casp3_F3_KO_stds = np.array([np.std(data) for data in DATA])

        ax[1,st+1].plot(time, Casp3_F3_KO_means, color=[0.0,0.8,0.0], lw=2, zorder=1)
        ax[1,st+1].scatter(time, Casp3_F3_KO_means, color=[0.0,0.6,0.0], edgecolor='k', s=50, zorder=2)
        ax[1,st+1].fill_between(time, Casp3_F3_KO_means-Casp3_F3_KO_stds, Casp3_F3_KO_means+Casp3_F3_KO_stds, alpha=0.1, color=[0.0,0.8,0.0])
        ax[1, st+1].sharey(ax[0, st+1])

        DATA = [np.array(Casp3_A12_WT[stage][n])/(np.array(A12_WT[n]) + np.array(Casp3_A12_WT[stage][n])) for n in range(len(A12_WT))]
        Casp3_A12_WT_means = np.array([np.mean(data) for data in DATA])
        Casp3_A12_WT_stds = np.array([np.std(data) for data in DATA])

        ax[0,st+1].plot(time, Casp3_A12_WT_means, color=[0.8,0.0,0.8], lw=2, zorder=1)
        ax[0,st+1].scatter(time, Casp3_A12_WT_means, color=[0.6,0.0,0.6], edgecolor='k', s=50, zorder=2)
        ax[0,st+1].fill_between(time, Casp3_A12_WT_means-Casp3_A12_WT_stds, Casp3_A12_WT_means+Casp3_A12_WT_stds, alpha=0.1, color=[0.8,0.0,0.8])

        ax[0, st+1].set_xticks([i for i in time])
        ax[0, st+1].set_xticklabels([None]*len(time))

        DATA = [np.array(Casp3_A12_KO[stage][n])/(np.array(A12_KO[n]) + np.array(Casp3_A12_KO[stage][n])) for n in range(len(A12_KO))]
        Casp3_A12_KO_means = np.array([np.mean(data) for data in DATA])
        Casp3_A12_KO_stds = np.array([np.std(data) for data in DATA])

        ax[1,st+1].plot(time, Casp3_A12_KO_means, color=[0.8,0.0,0.8], lw=2, zorder=1)
        ax[1,st+1].scatter(time, Casp3_A12_KO_means, color=[0.6,0.0,0.6], edgecolor='k', s=50, zorder=2)
        ax[1,st+1].fill_between(time, Casp3_A12_KO_means-Casp3_A12_KO_stds, Casp3_A12_KO_means+Casp3_A12_KO_stds, alpha=0.1, color=[0.8,0.0,0.8])
        ax[1, st+1].sharey(ax[0, st+1])
        
        ax[0, st+1].spines[['right', 'top']].set_visible(False)
        ax[1, st+1].spines[['right', 'top']].set_visible(False)
        ax[1, st+1].set_xticks([i for i in time])
        ax[1, st+1].set_xticklabels(TIMES)

    plt.tight_layout()
    plt.savefig(path_figure)
    plt.show()

    barwidth = 0.5

    fig, ax = plt.subplots(1, 4, figsize=(14, 3.5))

    ax[0].set_title("total cells")
    ax[1].set_title("early apoptotis")
    ax[2].set_title("mid apoptotis")
    ax[3].set_title("late apoptotis")
    ax[0].set_ylabel("cell number")
    ax[1].set_ylabel("apoptosis proportion")

    from scipy.stats import ttest_ind

    pvalues = []
    for cc, condition in enumerate(F3_WT):
        tt = ttest_ind(np.array(F3_WT[cc]), np.array(F3_KO[cc]))
        pvalues.append(tt.pvalue)
        print(tt.pvalue)

    data = [np.array([np.mean(data) for data in F3_WT]), np.array([np.mean(data) for data in F3_KO])]
    data_std = [np.array([np.std(data) for data in F3_WT]), np.array([np.std(data) for data in F3_KO])]
    maxz = np.max(np.array(data) + np.array(data_std))
    x = np.arange(len(TIMES)) * 1.5 + 1
    for i in range(2):
        if i==0:
            offset=-barwidth/2
            color = [0.0,0.3,0.0]
        elif i == 1:
            offset = +barwidth/2
            color = [0.0,0.8,0.0]

        rects = ax[0].bar(x + offset, np.array(data[i]), barwidth, label=CONDS[i], color=color, edgecolor="k")
        ax[0].errorbar(x+offset,  np.array(data[i]), yerr=np.array(data_std[i]), capsize=5, fmt="o", color="k")
        
        if i==1:
            for p, pval in enumerate(pvalues):
                x1 = x[p]-offset
                x2 = x[p]+offset
                y = maxz*1.1
                h = maxz*0.01
                ax[0].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c="k")
                if pval > 0.05:
                    ax[0].text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.01:
                    ax[0].text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.001:
                    ax[0].text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color="k", fontsize=12)
                else:
                    ax[0].text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color="k", fontsize=12)
                        
    ax[0].set_xticks(x, TIMES)
    ax[0].spines[['right', 'top']].set_visible(False)

    DATA1 = [np.array(Casp3_F3_WT[APO_STAGES[0]][n])/(np.array(F3_WT[n]) + np.array(Casp3_F3_WT[APO_STAGES[0]][n])) for n in range(len(F3_WT))]
    DATA2 = [np.array(Casp3_F3_KO[APO_STAGES[0]][n])/(np.array(F3_KO[n]) + np.array(Casp3_F3_KO[APO_STAGES[0]][n])) for n in range(len(F3_KO))]

    pvalues = []
    for cc, condition in enumerate(DATA1):
        tt = ttest_ind(np.array(DATA1[cc]), np.array(DATA2[cc]))
        pvalues.append(tt.pvalue)
        print(tt.pvalue)

    data = [np.array([np.mean(data) for data in DATA1]), np.array([np.mean(data) for data in DATA2])]
    data_std = [np.array([np.std(data) for data in DATA1]), np.array([np.std(data) for data in DATA2])]
    maxz = np.max(np.array(data) + np.array(data_std))
    x = np.arange(len(TIMES)) * 1.5 + 1
    for i in range(2):
        if i==0:
            offset=-barwidth/2
            color = [0.0,0.3,0.0]
        elif i == 1:
            offset = +barwidth/2
            color = [0.0,0.8,0.0]

        rects = ax[1].bar(x + offset, np.array(data[i]), barwidth, label=CONDS[i], color=color, edgecolor="k")
        ax[1].errorbar(x+offset,  np.array(data[i]), yerr=np.array(data_std[i]), capsize=5, fmt="o", color="k")
        
        if i==1:
            for p, pval in enumerate(pvalues):
                x1 = x[p]-offset
                x2 = x[p]+offset
                y = maxz*1.1
                h = maxz*0.01
                ax[1].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c="k")
                if pval > 0.05:
                    ax[1].text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.01:
                    ax[1].text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.001:
                    ax[1].text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color="k", fontsize=12)
                else:
                    ax[1].text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color="k", fontsize=12)
                        
    ax[1].set_xticks(x, TIMES)
    ax[1].spines[['right', 'top']].set_visible(False)
    
    DATA1 = [np.array(Casp3_F3_WT[APO_STAGES[1]][n])/(np.array(F3_WT[n]) + np.array(Casp3_F3_WT[APO_STAGES[1]][n])) for n in range(len(F3_WT))]
    DATA2 = [np.array(Casp3_F3_KO[APO_STAGES[1]][n])/(np.array(F3_KO[n]) + np.array(Casp3_F3_KO[APO_STAGES[1]][n])) for n in range(len(F3_KO))]

    pvalues = []
    for cc, condition in enumerate(DATA1):
        tt = ttest_ind(np.array(DATA1[cc]), np.array(DATA2[cc]))
        pvalues.append(tt.pvalue)
        print(tt.pvalue)

    data = [np.array([np.mean(data) for data in DATA1]), np.array([np.mean(data) for data in DATA2])]
    data_std = [np.array([np.std(data) for data in DATA1]), np.array([np.std(data) for data in DATA2])]
    maxz = np.max(np.array(data) + np.array(data_std))
    x = np.arange(len(TIMES)) * 1.5 + 1
    for i in range(2):
        if i==0:
            offset=-barwidth/2
            color = [0.0,0.3,0.0]
        elif i == 1:
            offset = +barwidth/2
            color = [0.0,0.8,0.0]

        rects = ax[2].bar(x + offset, np.array(data[i]), barwidth, label=CONDS[i], color=color, edgecolor="k")
        ax[2].errorbar(x+offset,  np.array(data[i]), yerr=np.array(data_std[i]), capsize=5, fmt="o", color="k")
        
        if i==1:
            for p, pval in enumerate(pvalues):
                x1 = x[p]-offset
                x2 = x[p]+offset
                y = maxz*1.1
                h = maxz*0.01
                ax[2].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c="k")
                if pval > 0.05:
                    ax[2].text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.01:
                    ax[2].text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.001:
                    ax[2].text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color="k", fontsize=12)
                else:
                    ax[2].text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color="k", fontsize=12)
                        
    ax[2].set_xticks(x, TIMES)
    ax[2].spines[['right', 'top']].set_visible(False)

    DATA1 = [np.array(Casp3_F3_WT[APO_STAGES[2]][n])/(np.array(F3_WT[n]) + np.array(Casp3_F3_WT[APO_STAGES[2]][n])) for n in range(len(F3_WT))]
    DATA2 = [np.array(Casp3_F3_KO[APO_STAGES[2]][n])/(np.array(F3_KO[n]) + np.array(Casp3_F3_KO[APO_STAGES[2]][n])) for n in range(len(F3_KO))]

    pvalues = []
    for cc, condition in enumerate(DATA1):
        tt = ttest_ind(np.array(DATA1[cc]), np.array(DATA2[cc]))
        pvalues.append(tt.pvalue)
        print(tt.pvalue)

    data = [np.array([np.mean(data) for data in DATA1]), np.array([np.mean(data) for data in DATA2])]
    data_std = [np.array([np.std(data) for data in DATA1]), np.array([np.std(data) for data in DATA2])]
    maxz = np.max(np.array(data) + np.array(data_std))
    x = np.arange(len(TIMES)) * 1.5 + 1
    for i in range(2):
        if i==0:
            offset=-barwidth/2
            color = [0.0,0.3,0.0]
        elif i == 1:
            offset = +barwidth/2
            color = [0.0,0.8,0.0]

        rects = ax[3].bar(x + offset, np.array(data[i]), barwidth, label=CONDS[i], color=color, edgecolor="k")
        ax[3].errorbar(x+offset,  np.array(data[i]), yerr=np.array(data_std[i]), capsize=5, fmt="o", color="k")
        
        if i==1:
            for p, pval in enumerate(pvalues):
                x1 = x[p]-offset
                x2 = x[p]+offset
                y = maxz*1.1
                h = maxz*0.01
                ax[3].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c="k")
                if pval > 0.05:
                    ax[3].text((x1+x2)*.5, y+h, "ns", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.01:
                    ax[3].text((x1+x2)*.5, y+h, "*", ha='center', va='bottom', color="k", fontsize=12)
                elif pval > 0.001:
                    ax[3].text((x1+x2)*.5, y+h, "**", ha='center', va='bottom', color="k", fontsize=12)
                else:
                    ax[3].text((x1+x2)*.5, y+h, "***", ha='center', va='bottom', color="k", fontsize=12)
                        
    ax[3].set_xticks(x, TIMES)
    ax[3].spines[['right', 'top']].set_visible(False)
    ax[3].legend()
    plt.tight_layout()
    path_figure = "{}counts_comparison_.pdf".format(path_figures_exp)
    plt.savefig(path_figure)
    
plt.show()



