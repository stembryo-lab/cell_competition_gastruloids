"""
Compute cell counts and apoptotic proportions from saved segmentations.

Loads the precomputed F3/A12 segmentations and Casp3 segmentations (early/mid/late),
counts objects per embryo, assigns apoptotic events to F3 vs A12 by channel intensity,
and exports per-file counts to a CSV for downstream plotting/statistics.
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

    all_files = []
    all_data = []
    all_names = ["file_name", "F3", "A12", "apo_F3_early", "apo_A12_early", "apo_F3_mid", "apo_A12_mid", "apo_F3_late", "apo_A12_late"]
    for T, TIME in enumerate(TIMES):
        for C,COND in enumerate(CONDS):
            path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/{}/'.format(EXP, TIME, COND)
            path_save_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/ctobjects/{}/{}/'.format(EXP, TIME, COND)
                
            ### GET FULL FILE NAME AND FILE CODE ###
            files = get_file_names(path_data_dir)

            all_files = [*all_files, *files]
            for f, file in enumerate(files):
                data = []
                if COND == "WT":
                    data.append(F3_WT[T][f])
                    data.append(A12_WT[T][f])
                    data.append(Casp3_F3_WT[APO_STAGES[0]][T][f])
                    data.append(Casp3_A12_WT[APO_STAGES[0]][T][f])
                    data.append(Casp3_F3_WT[APO_STAGES[1]][T][f])
                    data.append(Casp3_A12_WT[APO_STAGES[1]][T][f])
                    data.append(Casp3_F3_WT[APO_STAGES[2]][T][f])
                    data.append(Casp3_A12_WT[APO_STAGES[2]][T][f])
                elif COND == "KO":
                    data.append(F3_KO[T][f])
                    data.append(A12_KO[T][f])
                    data.append(Casp3_F3_KO[APO_STAGES[0]][T][f])
                    data.append(Casp3_A12_KO[APO_STAGES[0]][T][f])
                    data.append(Casp3_F3_KO[APO_STAGES[1]][T][f])
                    data.append(Casp3_A12_KO[APO_STAGES[1]][T][f])
                    data.append(Casp3_F3_KO[APO_STAGES[2]][T][f])
                    data.append(Casp3_A12_KO[APO_STAGES[2]][T][f])
                all_data.append(data)

    import csv
    # Output CSV file path
    output_file = "{}data_counts.csv".format(path_figures_exp)
    # Write to CSV
    with open(output_file, mode="w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(all_names)
        
        # Write the data rows
        for file, values in zip(all_files, all_data):
            csv_writer.writerow([file] + values)
