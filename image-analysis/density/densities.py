### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
TIMES = ["48hr", "72hr", "96hr"]
CONDS = ["WT", "KO"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

path_figures = "/home/pablo/Desktop/PhD/projects/GastruloidCompetition/results/densities/"

for EXP in EXPERIMENTS:
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)
    
    for number_of_neighs in [5 ,10, 15, 20, 30, 50, 75, 100, 200]:
        fig, ax = plt.subplots(2,3, figsize=(12,6), sharey=True, sharex='col')
        
        for ap, apo_stage in enumerate(["early", "mid", "late"]):

            densities_F3_all = [[[], [], []], [[], [], []]]
            densities_A12_all = [[[], [], []], [[], [], []]]
            densities_F3_all_apo = [[[], [], []], [[], [], []]]
            densities_A12_all_apo = [[[], [], []], [[], [], []]]
            gastruloid_sizes = [[[], [], []], [[], [], []]]
                
            for TTT, TIME in enumerate(TIMES):
                if EXP=="2023_11_17_Casp3":
                    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                    if TIME=="96hr":
                        channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
                elif EXP=="2024_03_Casp3":
                    channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                    
                for CCC, COND in enumerate(CONDS):
                    path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/{}/'.format(EXP, TIME, COND)
                    path_save_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/ctobjects/{}/{}/'.format(EXP, TIME, COND)
                    files = get_file_names(path_data_dir)

                    densities_F3 = []
                    densities_A12 = []
                    densities_F3_apo = []
                    densities_A12_apo = []
                    
                    for f, file in enumerate(files):
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

                        area = np.mean(areas)
                        dim = 2*np.sqrt(area/np.pi)
                                            
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
                            
                        Casp3_F3 = 0
                        Casp3_A12 = 0
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
                                Casp3_F3+=1
                            else:
                                fates.append(3)
                                Casp3_A12+=1
                            
                        centers = np.array(centers)
                        fates = np.array(fates)
                        
                        nbrs = NearestNeighbors(n_neighbors=number_of_neighs+1, algorithm='ball_tree').fit(centers)
                        distances, neighs = nbrs.kneighbors(centers)

                        dist_th = (dim*xyres)*1000.0 #microns
                        dist_th_near = (dim*xyres)*0.25
                        
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

                        densities = [1/(np.mean(dists))**3 for dists in true_dists]

                        _densities_F3 = [densities[n] for n in range(len(fates)) if fates[n] == 0]
                        _densities_A12 = [densities[n] for n in range(len(fates)) if fates[n] == 1]
                        _densities_F3_apo = [densities[n] for n in range(len(fates)) if fates[n] == 2]
                        _densities_A12_apo = [densities[n] for n in range(len(fates)) if fates[n] == 3]
                            
                        densities_F3 = [*densities_F3, *_densities_F3]
                        densities_A12 = [*densities_A12, *_densities_A12]
                        densities_F3_apo = [*densities_F3_apo, *_densities_F3_apo]
                        densities_A12_apo = [*densities_A12_apo, *_densities_A12_apo]

                    densities_F3_all[CCC][TTT] = [*densities_F3_all[CCC][TTT], *densities_F3]
                    densities_A12_all[CCC][TTT] = [*densities_A12_all[CCC][TTT], *densities_A12]
                    densities_F3_all_apo[CCC][TTT] = [*densities_F3_all_apo[CCC][TTT], *densities_F3_apo]
                    densities_A12_all_apo[CCC][TTT] = [*densities_A12_all_apo[CCC][TTT], *densities_A12_apo]
                    gastruloid_sizes[CCC][TTT] = [*gastruloid_sizes[CCC][TTT], len(centers)]

            F3s_mean_dens_WT = np.array([np.nanmean(densities_F3_all[0][i]) for i in range(3)])
            F3s_mean_dens_KO = np.array([np.nanmean(densities_F3_all[1][i]) for i in range(3)])
            F3s_mean_dens_apo_WT = np.array([np.nanmean(densities_F3_all_apo[0][i]) for i in range(3)])
            F3s_mean_dens_apo_KO = np.array([np.nanmean(densities_F3_all_apo[1][i]) for i in range(3)])

            A12s_mean_dens_WT = np.array([np.nanmean(densities_A12_all[0][i]) for i in range(3)])
            A12s_mean_dens_KO = np.array([np.nanmean(densities_A12_all[1][i]) for i in range(3)])
            A12s_mean_dens_apo_WT = np.array([np.nanmean(densities_A12_all_apo[0][i]) for i in range(3)])
            A12s_mean_dens_apo_KO = np.array([np.nanmean(densities_A12_all_apo[1][i]) for i in range(3)])

            g_sizes_WT = np.array([np.nanmean(gastruloid_sizes[0][i]) for i in range(3)])
            g_sizes_KO = np.array([np.nanmean(gastruloid_sizes[1][i]) for i in range(3)])

            F3s_stds_dens_WT = np.array([np.nanstd(densities_F3_all[0][i]) for i in range(3)])
            F3s_stds_dens_KO = np.array([np.nanstd(densities_F3_all[1][i]) for i in range(3)])
            A12s_stds_dens_WT = np.array([np.nanstd(densities_A12_all[0][i]) for i in range(3)])
            A12s_stds_dens_KO = np.array([np.nanstd(densities_A12_all[1][i]) for i in range(3)])
            g_sizes_stds_WT = np.array([np.nanstd(gastruloid_sizes[0][i]) for i in range(3)])
            g_sizes_stds_KO = np.array([np.nanstd(gastruloid_sizes[1][i]) for i in range(3)])

            conds = np.array([48, 72, 96])

            ax[0, ap].set_title("{} apoptotis".format(apo_stage))
            ax[0, ap].plot(conds, F3s_mean_dens_WT, color="green", lw=3, label="F3")
            ax[0, ap].scatter(conds, F3s_mean_dens_WT, color="green", s=100, edgecolor="k", zorder=10)
            ax[0, ap].plot(conds, F3s_mean_dens_apo_WT, color="green", lw=3, ls='--', label="F3 - apo")
            ax[0, ap].scatter(conds, F3s_mean_dens_apo_WT, color="green", s=100, edgecolor="k", zorder=10)

            ax[0, ap].plot(conds, A12s_mean_dens_WT, color="magenta", lw=3, label="A12")
            ax[0, ap].scatter(conds, A12s_mean_dens_WT, color="magenta", s=100, edgecolor="k", zorder=10)
            ax[0, ap].plot(conds, A12s_mean_dens_apo_WT, color="magenta", lw=3, ls='--', label="A12 - apo")
            ax[0, ap].scatter(conds, A12s_mean_dens_apo_WT, color="magenta", s=100, edgecolor="k", zorder=10)

            ax[1, ap].plot(conds, F3s_mean_dens_KO, color="green", lw=3, label="F3")
            ax[1, ap].scatter(conds, F3s_mean_dens_KO, color="green", s=100, edgecolor="k", zorder=10)
            ax[1, ap].plot(conds, F3s_mean_dens_apo_KO, color="green", lw=3, ls='--', label="F3 - apo")
            ax[1, ap].scatter(conds, F3s_mean_dens_apo_KO, color="green", s=100, edgecolor="k", zorder=10)

            ax[1, ap].plot(conds, A12s_mean_dens_KO, color="magenta", lw=3, label="A12")
            ax[1, ap].scatter(conds, A12s_mean_dens_KO, color="magenta", s=100, edgecolor="k", zorder=10)
            ax[1, ap].plot(conds, A12s_mean_dens_apo_KO, color="magenta", lw=3, ls='--', label="A12 - apo")
            ax[1, ap].scatter(conds, A12s_mean_dens_apo_KO, color="magenta", s=100, edgecolor="k", zorder=10)

            ax[1, ap].set_xticks(conds)
            ax[1, ap].set_xlabel("Time (hr)")
            if ap==0:
                ax[0, ap].set_ylabel(r"mean $\rho_\mathrm{{local}}$  $\mu \mathrm{{m}}^{{-3}}$")
                ax[1, ap].set_ylabel(r"mean $\rho_\mathrm{{local}}$  $\mu \mathrm{{m}}^{{-3}}$")

            ax[0, ap].spines[['right', 'top']].set_visible(False)
            ax[1, ap].spines[['right', 'top']].set_visible(False)

            miny = np.min([*F3s_mean_dens_KO,*A12s_mean_dens_KO,*F3s_mean_dens_apo_KO, *A12s_mean_dens_apo_KO, *F3s_mean_dens_WT,*A12s_mean_dens_WT])
            maxy = np.max([*F3s_mean_dens_KO,*A12s_mean_dens_KO,*F3s_mean_dens_apo_KO, *A12s_mean_dens_apo_KO, *F3s_mean_dens_WT,*A12s_mean_dens_WT])
            
        plt.tight_layout()
        plt.savefig(path_figures_exp+"density_{}.svg".format(number_of_neighs))
        plt.savefig(path_figures_exp+"density_{}.pdf".format(number_of_neighs))

    plt.show()