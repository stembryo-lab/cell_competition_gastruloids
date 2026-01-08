"""
Detect and remove misclassified nuclei in the p53-degron dataset.

For each z-plane, quantify lineage-marker intensities (F3/mCherry and A12/emiRFP)
inside nuclear masks, cluster cells in log-intensity space (k=2), and flag
“spillover” assignments (cluster ≠ segmentation label and/or low-signal tail).
Then reload the saved segmentations and delete the flagged labels, updating
labels in-place. Also saves per-z diagnostic scatter plots.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels


import numpy as np
import matplotlib.pyplot as plt

F3_all = [[] for z in range(10)]
F3_F3 = [[] for z in range(10)]
F3_A12 = [[] for z in range(10)]
F3_DAPI = [[] for z in range(10)]

A12_all = [[] for z in range(10)]
A12_F3 = [[] for z in range(10)]
A12_A12 = [[] for z in range(10)]
A12_DAPI = [[] for z in range(10)]

DAPI_all = [[] for z in range(10)]

colors = [[] for z in range(10)]
fates = [[] for z in range(10)]

CONDS = ["auxin_48-72_48", "auxin_48-72_72" , "auxin_48-72_96", "auxin_72-96_72", "auxin_72-96_96", "noauxin_72", "noauxin_96", "secondaryonly"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

zs = []
all_files = []

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53_degron/segmentation_results/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53_degron/figures/"
check_or_create_dir(path_to_save_results)

for COND in CONDS:
    path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/2025_09_09_OsTIRMosaic_p53Timecourse/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    files = get_file_names(path_data_dir)

    channel_names = ["A12", "p53", "F3", "DAPI"]

    for f, file in enumerate(files):
                
        all_files.append(file)
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

        zs.append(CT_A12.hyperstack.shape[1])
        
        ch_F3 = channel_names.index("F3")
        ch_A12 = channel_names.index("A12")
        ch_DAPI = channel_names.index("DAPI")

        for cell in CT_F3.jitcells:
            center = cell.centers[0]
            z = int(center[0])
            zid = cell.zs[0].index(z)
            mask = cell.masks[0][zid]

            F3_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_F3,:,:][mask[:,1], mask[:,0]]))
            A12_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_A12,:,:][mask[:,1], mask[:,0]]))
            DAPI_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_DAPI,:,:][mask[:,1], mask[:,0]]))

            F3_F3[z].append(np.mean(CT_A12.hyperstack[0,z,ch_F3,:,:][mask[:,1], mask[:,0]]))
            F3_A12[z].append(np.mean(CT_A12.hyperstack[0,z,ch_A12,:,:][mask[:,1], mask[:,0]]))
            F3_DAPI[z].append(np.mean(CT_A12.hyperstack[0,z,ch_DAPI,:,:][mask[:,1], mask[:,0]]))

            colors[z].append([0.0,0.8,0.0, 0.3])
            fates[z].append("F3")

        for cell in CT_A12.jitcells:
            center = cell.centers[0]
            z = int(center[0])
            zid = cell.zs[0].index(z)
            mask = cell.masks[0][zid]

            F3_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_F3,:,:][mask[:,1], mask[:,0]]))
            A12_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_A12,:,:][mask[:,1], mask[:,0]]))
            DAPI_all[z].append(np.mean(CT_A12.hyperstack[0,z,ch_DAPI,:,:][mask[:,1], mask[:,0]]))
            
            A12_F3[z].append(np.mean(CT_A12.hyperstack[0,z,ch_F3,:,:][mask[:,1], mask[:,0]]))
            A12_A12[z].append(np.mean(CT_A12.hyperstack[0,z,ch_A12,:,:][mask[:,1], mask[:,0]]))
            A12_DAPI[z].append(np.mean(CT_A12.hyperstack[0,z,ch_DAPI,:,:][mask[:,1], mask[:,0]]))                
            colors[z].append([0.8,0.0,0.8, 0.3])
            fates[z].append("A12")

F3_means = np.array([np.mean(f3) for f3 in F3_F3])
F3_stds = np.array([np.std(f3) for f3 in F3_F3])

A12_means = np.array([np.mean(a12) for a12 in A12_A12])
A12_stds = np.array([np.std(a12) for a12 in A12_A12])

DAPI_means = np.array([np.mean(dapi) for dapi in DAPI_all])
DAPI_stds = np.array([np.std(dapi) for dapi in DAPI_all])

fig, ax = plt.subplots()

ax.plot(range(10), F3_means, color=[0.0,0.8,0.0, 1.0], label="H2B-mCherry on F3")
ax.fill_between(range(10), F3_means - F3_stds, F3_means + F3_stds, color=[0.0,0.8,0.0], alpha=0.2)

ax.plot(range(10), A12_means, color=[0.8,0.0,0.8, 1.0], label="H2B-emiRFP on A12")
ax.fill_between(range(10), A12_means - A12_stds, A12_means + A12_stds, color=[0.8,0.0,0.8], alpha=0.2)

ax.plot(range(10), DAPI_means, color="cyan", label="DAPI on all cells")
ax.fill_between(range(10), DAPI_means - DAPI_stds, DAPI_means + DAPI_stds, color="cyan", alpha=0.2)

ax.set_xlabel("z")
ax.set_ylabel("fluoro [a.u.]")
ax.title("Mean with std ribbon")
ax.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

remove_cell = [[] for z in range(10)]

for z in range(10):
    data1 = F3_all[z]
    data2 = A12_all[z]
    checkvar = fates[z]
    X = np.transpose(np.asarray([np.log(data1), np.log(data2)]))

    kmeans = KMeans(n_clusters=2)
    kmeans.fit(X)
    cluster_centers = kmeans.cluster_centers_
    clustered_fates = np.argsort(cluster_centers[:,0])
    labels = kmeans.labels_

    colors_clustering = []

    A12th = np.percentile(data2,20)
    F3th = np.percentile(data1,20)

    for i, lab in enumerate(labels):
        if lab==clustered_fates[0]:
            if checkvar[i]!="A12":
                colors_clustering.append([0.2, 0.2, 0.2, 1.0])
                remove_cell[z].append(True)
            else:
                colors_clustering.append([0.8, 0.0, 0.8, 0.7])
                remove_cell[z].append(False)
        elif lab==clustered_fates[-1]:
            if checkvar[i]!="F3":
                colors_clustering.append([0.2, 0.2, 0.2, 1.0])
                remove_cell[z].append(True)
            else:
                remove_cell[z].append(False)
                colors_clustering.append([0.0, 0.8, 0.0, 0.7])
        else:
            remove_cell[z].append(False)
            colors_clustering.append([0.0, 0.0, 0.0, 0.7])

        if checkvar[i]=="A12":
            if data2[i] < A12th:
                remove_cell[z][-1]=True
                colors_clustering[-1] = [0.2, 0.2, 0.2, 1.0]
        else:
            if data1[i] < F3th:
                remove_cell[z][-1]=True
                colors_clustering[-1] = [0.2, 0.2, 0.2, 1.0]
            
    # Plot the original data points and cluster centers
    from matplotlib.lines import Line2D
    plt.figure(figsize=(8, 6))
    plt.scatter(X[:, 0], X[:, 1], c=colors_clustering, edgecolors='k')
    plt.xlabel('log(emiRFG)')
    plt.ylabel('log(mCherry)')
    plt.tight_layout()

    # Custom legend handles
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="A12",
            markerfacecolor=(0.8, 0.0, 0.8, 0.7), markersize=10),
        Line2D([0], [0], marker="o", color="w", label="F3",
            markerfacecolor=(0.0, 0.8, 0.0, 0.7), markersize=10),
        Line2D([0], [0], marker="o", color="w", label="spillover cells",
            markerfacecolor=(0.2, 0.2, 0.2, 1.0), markersize=10),
    ]
    plt.legend(handles=legend_elements, loc="best", frameon=True)
    plt.savefig("{}clusteringloglog_z{}.svg".format(path_to_save_results, z))
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.scatter(data1, data2, c=colors_clustering, edgecolors='k')
    plt.xlabel('emiRFP')
    plt.ylabel('mCherry')
    plt.tight_layout()

    # Custom legend handles
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="A12",
            markerfacecolor=(0.8, 0.0, 0.8, 0.7), markersize=10),
        Line2D([0], [0], marker="o", color="w", label="F3",
            markerfacecolor=(0.0, 0.8, 0.0, 0.7), markersize=10),
        Line2D([0], [0], marker="o", color="w", label="spillover cells",
            markerfacecolor=(0.2, 0.2, 0.2, 1.0), markersize=10),
    ]
    plt.legend(handles=legend_elements, loc="best", frameon=True)
    plt.savefig("{}clustering_z{}.svg".format(path_to_save_results, z))

plt.show()

# Now is time for the removal
current_zid = [0 for z in range(10)]
for COND in CONDS:
    path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/2025_09_09_OsTIRMosaic_p53Timecourse/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    files = get_file_names(path_data_dir)

    channel_names = ["A12", "p53", "F3", "DAPI"]

    for f, file in enumerate(files):
                
        all_files.append(file)
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

        labs_to_rem = []
        for cell in CT_F3.jitcells:
            center = cell.centers[0]
            z = int(center[0])
            if remove_cell[z][current_zid[z]]:
                labs_to_rem.append(cell.label)
            current_zid[z]+=1
        
        for lab in labs_to_rem:
            CT_F3._del_cell(lab)  
            
        CT_F3.update_labels() 
    
        labs_to_rem = []
        for cell in CT_A12.jitcells:
            center = cell.centers[0]
            z = int(center[0])
            if remove_cell[z][current_zid[z]]:
                labs_to_rem.append(cell.label)
            current_zid[z]+=1
        
        for lab in labs_to_rem:
            CT_A12._del_cell(lab)  
            
        CT_A12.update_labels() 
