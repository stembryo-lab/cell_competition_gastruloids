### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

### LOAD STARDIST MODEL ###
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained('2D_versatile_fluo')

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

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
TIMES = ["48hr", "72hr", "96hr"]
CONDS = ["WT", "KO"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/neighborhood/"
check_or_create_dir(path_figures)

for EXP in EXPERIMENTS:
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)
    
    for number_of_neighs in [5 ,10, 15, 20, 30, 50, 75, 100, 200]:
        # check_or_create_dir(path_figures+"{}/".format(number_of_neighs))
        fig, ax = plt.subplots(2,3, figsize=(12,6), sharey=True, sharex='col')
        for ap, apo_stage in enumerate(["early", "mid", "late"]):
            
            # check_or_create_dir(path_figures+"{}/{}/".format(number_of_neighs, apo_stage))

            NEIGHS_F3 = [[], []]
            NEIGHS_A12 = [[], []]
            NEIGHS_F3_apo = [[], []]
            NEIGHS_A12_apo = [[], []]
                
            for TTT, TIME in enumerate(TIMES):
                if EXP=="2023_11_17_Casp3":
                    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                    if TIME=="96hr":
                        channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
                elif EXP=="2024_03_Casp3":
                    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                    
                for CCC, COND in enumerate(CONDS):
                    path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
                    path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
                    files_data = get_file_names(path_data_dir)

                    neighs_fates_F3_sum = np.zeros((len(files_data), 2))
                    neighs_fates_A12_sum = np.zeros((len(files_data), 2))
                    neighs_fates_Casp3_F3_sum = np.zeros((len(files_data), 2))
                    neighs_fates_Casp3_A12_sum = np.zeros((len(files_data), 2))

                    total_f3s = []
                    total_a12s = []
                
                    for f, file in enumerate(files_data):
                        path_data = path_data_dir+file
                        file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
                        path_save = correct_path(path_save_dir+embcode)
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


                        ### CORRECT DISTRIBUTIONS ###

                        F3_dist = []
                        areas = []
                        for cell in CT_F3.jitcells: 
                            for zid, z in enumerate(cell.zs[0]):
                                mask = cell.masks[0][zid]
                                img = CT_F3.hyperstack[0,z, channel_names.index("F3")]
                                F3_dist.append(np.mean(img[mask[:,1], mask[:,0]]))

                        A12_dist = []
                        for cell in CT_A12.jitcells: 
                            for zid, z in enumerate(cell.zs[0]):
                                mask = cell.masks[0][zid]
                                areas.append(len(mask))
                                img = CT_A12.hyperstack[0,z, channel_names.index("A12")]
                                A12_dist.append(np.mean(img[mask[:,1], mask[:,0]]))

                        total_f3s.append(len(CT_F3.jitcells))
                        total_a12s.append(len(CT_A12.jitcells))

                        area = np.mean(areas)
                        dim = 2*np.sqrt(area/np.pi)
                                        
                        ## Now contract the shape as much as we want. 
                        F3_dist = np.array(F3_dist)
                        A12_dist = np.array(A12_dist)
                        
                        mdiff = np.mean(F3_dist) - np.mean(A12_dist)
                        
                        zres = CT_F3.metadata["Zresolution"]
                        xyres = CT_F3.metadata["XYresolution"]
                        
                        fates = []
                        centers = []
                        for cell in CT_F3.jitcells:
                            fates.append(0)
                            centers.append(cell.centers[0]*[zres, xyres, xyres])
                        for cell in CT_A12.jitcells:
                            fates.append(1)
                            centers.append(cell.centers[0]*[zres, xyres, xyres])

                        len_pre_casp3 = len(centers)
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

                            centers.append(cell.centers[0]*[zres, xyres, xyres])

                            zz = np.int64(cell.centers[0][0])
                            idx = cell.zs[0].index(zz)
                            
                            cell_f3 = f3[idx]
                            cell_a12 = a12[idx]

                            if mdiff > 0:
                                cell_a12 += mdiff
                            else: 
                                cell_f3 -= mdiff 
                                
                            if cell_f3 > cell_a12:
                                fates.append(2)
                            else:
                                fates.append(3)
                            
                        centers = np.array(centers)
                        fates = np.array(fates)
                        
                        nbrs = NearestNeighbors(n_neighbors=number_of_neighs+1, algorithm='ball_tree').fit(centers)
                        distances, neighs = nbrs.kneighbors(centers)

                        dist_th = (dim*xyres)*10000.0 #microns
                        dist_th_near = (dim*xyres)
                        
                        true_neighs = []
                        true_dists = []
                        for p, neigh_p in enumerate(neighs):
                            true_neigh_p = []
                            true_dists_p = []
                            for neigh in neigh_p[1:]:
                                dist = np.linalg.norm(centers[p]-centers[neigh])
                                if dist < dist_th:
                                    if dist > dist_th_near:
                                        if (neigh < len_pre_casp3) or (apo_stage=="early"):
                                            true_dists_p.append(dist)
                                            true_neigh_p.append(neigh)
                                if len(true_neigh_p) == number_of_neighs: break
                            true_dists.append(true_dists_p)
                            true_neighs.append(true_neigh_p)
                        
                        neighs_fates = []
                        for p, neigh_p in enumerate(true_neighs):
                            lbs = []
                            fts = []
                            for neigh in neigh_p:
                                fts.append(fates[neigh])
                            neighs_fates.append(fts)

                        neighs_fates_F3 = [n for i,n in enumerate(neighs_fates) if fates[i] == 0]
                        neighs_fates_A12 = [n for i,n in enumerate(neighs_fates) if fates[i] == 1]
                        neighs_fates_Casp3_F3 = [n for i,n in enumerate(neighs_fates) if fates[i] == 2]
                        neighs_fates_Casp3_A12 = [n for i,n in enumerate(neighs_fates) if fates[i] == 3]
                        
                        for n_fates in neighs_fates_F3:
                            for _f in n_fates:
                                if _f in [0]:
                                    neighs_fates_F3_sum[f,0] += 1
                                elif _f in [1]:
                                    neighs_fates_F3_sum[f,1] += 1
                                    
                        for n_fates in neighs_fates_A12:
                            for _f in n_fates:
                                if _f in [0]:
                                    neighs_fates_A12_sum[f,0] += 1
                                elif _f in [1]:
                                    neighs_fates_A12_sum[f,1] += 1

                        for n_fates in neighs_fates_Casp3_F3:
                            for _f in n_fates:
                                if _f in [0]:
                                    neighs_fates_Casp3_F3_sum[f,0] += 1
                                elif _f in [1]:
                                    neighs_fates_Casp3_F3_sum[f,1] += 1

                        for n_fates in neighs_fates_Casp3_A12:
                            for _f in n_fates:
                                if _f in [0]:
                                    neighs_fates_Casp3_A12_sum[f,0] += 1
                                elif _f in [1]:
                                    neighs_fates_Casp3_A12_sum[f,1] += 1

                        f3_a12_norm = np.array([total_f3s[f], total_a12s[f]])
                        
                        neighs_fates_A12_sum[f] /= f3_a12_norm
                        neighs_fates_A12_sum[f] /= np.sum(neighs_fates_A12_sum[f])
                        
                        neighs_fates_F3_sum[f] /= f3_a12_norm
                        neighs_fates_F3_sum[f] /= np.sum(neighs_fates_F3_sum[f])

                        neighs_fates_Casp3_F3_sum[f] /= f3_a12_norm
                        neighs_fates_Casp3_F3_sum[f] /= np.sum(neighs_fates_Casp3_F3_sum[f])

                        neighs_fates_Casp3_A12_sum[f] /= f3_a12_norm
                        neighs_fates_Casp3_A12_sum[f] /= np.sum(neighs_fates_Casp3_A12_sum[f])
                        
                    NEIGHS_F3[CCC].append(neighs_fates_F3_sum)
                    NEIGHS_F3_apo[CCC].append(neighs_fates_Casp3_F3_sum)
                    NEIGHS_A12[CCC].append(neighs_fates_A12_sum)
                    NEIGHS_A12_apo[CCC].append(neighs_fates_Casp3_A12_sum)
            
            F3s_mean_F3_neigh_WT = np.array([np.nanmean(NEIGHS_F3[0][i][:, 0]) for i in range(len(TIMES))])
            F3s_mean_F3_neigh_KO = np.array([np.nanmean(NEIGHS_F3[1][i][:, 0]) for i in range(len(TIMES))])
            F3s_mean_A12_neigh_WT = np.array([np.nanmean(NEIGHS_F3[0][i][:, 1]) for i in range(len(TIMES))])
            F3s_mean_A12_neigh_KO = np.array([np.nanmean(NEIGHS_F3[1][i][:, 1]) for i in range(len(TIMES))])

            F3s_apo_mean_F3_neigh_WT = np.array([np.nanmean(NEIGHS_F3_apo[0][i][:, 0]) for i in range(len(TIMES))])
            F3s_apo_mean_F3_neigh_KO = np.array([np.nanmean(NEIGHS_F3_apo[1][i][:, 0]) for i in range(len(TIMES))])
            F3s_apo_mean_A12_neigh_WT = np.array([np.nanmean(NEIGHS_F3_apo[0][i][:, 1]) for i in range(len(TIMES))])
            F3s_apo_mean_A12_neigh_KO = np.array([np.nanmean(NEIGHS_F3_apo[1][i][:, 1]) for i in range(len(TIMES))])

            A12s_mean_F3_neigh_WT = np.array([np.nanmean(NEIGHS_A12[0][i][:, 0]) for i in range(len(TIMES))])
            A12s_mean_F3_neigh_KO = np.array([np.nanmean(NEIGHS_A12[1][i][:, 0]) for i in range(len(TIMES))])
            A12s_mean_A12_neigh_WT = np.array([np.nanmean(NEIGHS_A12[0][i][:, 1]) for i in range(len(TIMES))])
            A12s_mean_A12_neigh_KO = np.array([np.nanmean(NEIGHS_A12[1][i][:, 1]) for i in range(len(TIMES))])

            A12s_apo_mean_F3_neigh_WT = np.array([np.nanmean(NEIGHS_A12_apo[0][i][:, 0]) for i in range(len(TIMES))])
            A12s_apo_mean_F3_neigh_KO = np.array([np.nanmean(NEIGHS_A12_apo[1][i][:, 0]) for i in range(len(TIMES))])
            A12s_apo_mean_A12_neigh_WT = np.array([np.nanmean(NEIGHS_A12_apo[0][i][:, 1]) for i in range(len(TIMES))])
            A12s_apo_mean_A12_neigh_KO = np.array([np.nanmean(NEIGHS_A12_apo[1][i][:, 1]) for i in range(len(TIMES))])


            F3s_std_F3_neigh_WT = np.array([np.nanstd(NEIGHS_F3[0][i][:, 0]) for i in range(len(TIMES))])
            F3s_std_F3_neigh_KO = np.array([np.nanstd(NEIGHS_F3[1][i][:, 0]) for i in range(len(TIMES))])
            F3s_std_A12_neigh_WT = np.array([np.nanstd(NEIGHS_F3[0][i][:, 1]) for i in range(len(TIMES))])
            F3s_std_A12_neigh_KO = np.array([np.nanstd(NEIGHS_F3[1][i][:, 1]) for i in range(len(TIMES))])

            F3s_apo_std_F3_neigh_WT = np.array([np.nanstd(NEIGHS_F3_apo[0][i][:, 0]) for i in range(len(TIMES))])
            F3s_apo_std_F3_neigh_KO = np.array([np.nanstd(NEIGHS_F3_apo[1][i][:, 0]) for i in range(len(TIMES))])
            F3s_apo_std_A12_neigh_WT = np.array([np.nanstd(NEIGHS_F3_apo[0][i][:, 1]) for i in range(len(TIMES))])
            F3s_apo_std_A12_neigh_KO = np.array([np.nanstd(NEIGHS_F3_apo[1][i][:, 1]) for i in range(len(TIMES))])

            A12s_std_F3_neigh_WT = np.array([np.nanstd(NEIGHS_A12[0][i][:, 0]) for i in range(len(TIMES))])
            A12s_std_F3_neigh_KO = np.array([np.nanstd(NEIGHS_A12[1][i][:, 0]) for i in range(len(TIMES))])
            A12s_std_A12_neigh_WT = np.array([np.nanstd(NEIGHS_A12[0][i][:, 1]) for i in range(len(TIMES))])
            A12s_std_A12_neigh_KO = np.array([np.nanstd(NEIGHS_A12[1][i][:, 1]) for i in range(len(TIMES))])

            A12s_apo_std_F3_neigh_WT = np.array([np.nanstd(NEIGHS_A12_apo[0][i][:, 0]) for i in range(len(TIMES))])
            A12s_apo_std_F3_neigh_KO = np.array([np.nanstd(NEIGHS_A12_apo[1][i][:, 0]) for i in range(len(TIMES))])
            A12s_apo_std_A12_neigh_WT = np.array([np.nanstd(NEIGHS_A12_apo[0][i][:, 1]) for i in range(len(TIMES))])
            A12s_apo_std_A12_neigh_KO = np.array([np.nanstd(NEIGHS_A12_apo[1][i][:, 1]) for i in range(len(TIMES))])

            conds = np.array([48, 72, 96])

            ax[0, ap].set_title("{} apoptotis".format(apo_stage))
            ax[0, ap].plot(conds, F3s_mean_A12_neigh_WT, color=[0.9,0.0,0.9], lw=4, label="F3")
            ax[0, ap].scatter(conds, F3s_mean_A12_neigh_WT, color=[0.9,0.0,0.9], s=100, edgecolor="k", zorder=10)
            ax[0, ap].plot(conds, F3s_apo_mean_A12_neigh_WT, color=[0.5,0.0,0.5], lw=4, ls='--', label="F3 - apo")
            ax[0, ap].scatter(conds, F3s_apo_mean_A12_neigh_WT, color=[0.5,0.0,0.5], s=100, edgecolor="k", zorder=10)

            ax[1, ap].plot(conds, F3s_mean_A12_neigh_KO, color=[0.9,0.0,0.9], lw=4, label="F3")
            ax[1, ap].scatter(conds, F3s_mean_A12_neigh_KO, color=[0.9,0.0,0.9], s=100, edgecolor="k", zorder=10)
            ax[1, ap].plot(conds, F3s_apo_mean_A12_neigh_KO, color=[0.5,0.0,0.5], lw=4, ls='--', label="F3 - apo")
            ax[1, ap].scatter(conds, F3s_apo_mean_A12_neigh_KO, color=[0.5,0.0,0.5], s=100, edgecolor="k", zorder=10)

            ax[1, ap].set_xticks(conds)
            ax[1, ap].set_xlabel("Time (hr)")
            if ap==0:
                ax[0, ap].set_ylabel(r"$\overline{F_{A12-WT}}$ around F3s")
                ax[1, ap].set_ylabel(r"$\overline{F_{A12-KO}}$ around F3s")

            ax[0, ap].spines[['right', 'top']].set_visible(False)
            ax[1, ap].spines[['right', 'top']].set_visible(False)
        
            if ap==2:
                ax[1, ap].legend()
            
        plt.tight_layout()
        plt.savefig(path_figures_exp+"neighborhood_{}.svg".format(number_of_neighs))
        plt.savefig(path_figures_exp+"neighborhood_{}.pdf".format(number_of_neighs))
