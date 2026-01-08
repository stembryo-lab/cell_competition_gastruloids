"""
Compute local cell densities using k-nearest neighbors and export to CSV.

Loads saved segmentations (F3, A12, and Casp3 early/mid/late), builds a 3D kNN
graph from cell centers, estimates per-cell density from neighbor distances, and
writes per-embryo mean densities (by population/apoptotic stage) for multiple k.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels
import numpy as np
import csv
from sklearn.neighbors import NearestNeighbors

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
TIMES = ["48hr", "72hr", "96hr"]
CONDS = ["WT", "KO"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
path_figures = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/densities/"
check_or_create_dir(path_figures)

for EXP in EXPERIMENTS:
    path_figures_exp = path_figures+"{}/".format(EXP)
    check_or_create_dir(path_figures_exp)
    
    for number_of_neighs in [5 ,10, 15, 20, 30, 50, 75, 100, 200]:
        filenames_early = []
        filenames_mid = []
        filenames_late = []

        densities_F3 = []
        densities_A12 = []
        densities_F3_early_apo = []
        densities_A12_early_apo = []
        densities_F3_mid_apo = []
        densities_A12_mid_apo = []
        densities_F3_late_apo = []
        densities_A12_late_apo = []
        
        for ap, apo_stage in enumerate(["early", "mid", "late"]):
                
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
                    files = get_file_names(path_data_dir)

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
                        
                        if apo_stage=="early":
                            filenames_early.append(file)
                            densities_F3.append(np.nanmean(_densities_F3))
                            densities_A12.append(np.nanmean(_densities_A12))
                            densities_F3_early_apo.append(np.nanmean(_densities_F3_apo))
                            densities_A12_early_apo.append(np.nanmean(_densities_A12_apo))
                        elif apo_stage=="mid":
                            filenames_mid.append(file)
                            densities_F3_mid_apo.append(np.nanmean(_densities_F3_apo))
                            densities_A12_mid_apo.append(np.nanmean(_densities_A12_apo))
                        elif apo_stage=="late":
                            filenames_late.append(file)
                            densities_F3_late_apo.append(np.nanmean(_densities_F3_apo))
                            densities_A12_late_apo.append(np.nanmean(_densities_A12_apo))

        output_file = path_figures_exp+"densities_{}neighs.csv".format(number_of_neighs)

        colnames = ["file", "F3", "A12", "apo_F3_early", "apo_A12_early", "apo_F3_mid", "apo_A12_mid", "apo_F3_late", "apo_A12_late"]
        full_data = [colnames]
        for f in range(len(filenames_early)):
            dat = [filenames_early[f], densities_F3[f], densities_A12[f], densities_F3_early_apo[f], densities_A12_early_apo[f], densities_F3_mid_apo[f], densities_A12_mid_apo[f], densities_F3_late_apo[f], densities_A12_late_apo[f]]
            full_data.append(dat)

        # Output CSV file path

        # Write to CSV
        with open(output_file, mode="w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header
            csv_writer.writerows(full_data)
