"""
Compute radial cell distribution metrics from saved segmentations.

Loads F3/A12/Casp3 segmentations, estimates the gastruloid centroid, and segments
the gastruloid boundary using active contours (Chan-Vese / “active contours without edges”).
For each cell, computes distances to the centroid and to the closest boundary point,
and saves per-embryo distance arrays for downstream radial-distribution analysis.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels, compute_distance_xyz_jit, compute_dists_jit, EmbryoSegmentation, tif_reader_5D
import numpy as np

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
TIMES = ["48hr", "72hr", "96hr"]
CONDS = ["WT", "KO"]
apo_stages = ["early", "mid", "late"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
path_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/radial_distribution/values/"
check_or_create_dir(path_results)

for EXP in EXPERIMENTS:
    path_results_exp = path_results+"{}/".format(EXP)
    check_or_create_dir(path_results_exp)

    for TIME in TIMES:
        print(TIME)
        if EXP=="2023_11_17_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
            if TIME=="96hr":
                channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        elif EXP=="2024_03_Casp3":
            channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

        for COND in CONDS:
            print(COND)
            if EXP=="2023_11_17_Casp3":
                if TIME=="96hr":
                    if COND=="KO":
                        binths = [6,10,10,10]
                    elif COND=="WT":
                        binths = [12,10,5,6]                    
            elif EXP=="2024_03_Casp3":
                if TIME=="96hr":
                    if COND=="KO":
                        binths = [[5,4],[5.0, 4.0],[6.0, 5.0]]
                    elif COND=="WT":
                        binths = [[7.5, 5.0],[7.5, 5.0],[7.5, 5.0]]
            
            path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
            path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
            files_data = get_file_names(path_data_dir)
    
            ch_F3 = channel_names.index("F3")
            ch_A12 = channel_names.index("A12")
            ch_Casp3 = channel_names.index("Casp3")
            
            for apo_stage in apo_stages:
                path_save_results='{}{}_apoptosis/{}/{}/'.format(path_results_exp, apo_stage, TIME, COND)
                check_or_create_dir(path_save_results)

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

                    import numpy as np

                    ### CORRECT DISTRIBUTIONS ###

                    areas = []

                    F3_dist = []
                    for cell in CT_F3.jitcells: 
                        for zid, z in enumerate(cell.zs[0]):
                            mask = cell.masks[0][zid]
                            img = CT_F3.hyperstack[0,z, channel_names.index("F3")]
                            F3_dist.append(np.mean(img[mask[:,1], mask[:,0]]))
                            areas.append(len(mask))
                            
                    A12_dist = []
                    for cell in CT_A12.jitcells: 
                        for zid, z in enumerate(cell.zs[0]):
                            mask = cell.masks[0][zid]
                            areas.append(len(mask))
                            img = CT_A12.hyperstack[0,z, channel_names.index("A12")]
                            A12_dist.append(np.mean(img[mask[:,1], mask[:,0]]))

                    area = np.mean(areas)
                    dim = 2*np.sqrt(area/np.pi)

                    ## Now contract the shape as much as we want. 
                    F3_dist = np.array(F3_dist)
                    A12_dist = np.array(A12_dist)

                    mdiff = np.mean(F3_dist) - np.mean(A12_dist)
                    if mdiff > 0:
                        A12_dist += mdiff
                    else: 
                        F3_dist -= mdiff 

                    zres = CT_F3.metadata["Zresolution"]
                    xyres = CT_F3.metadata["XYresolution"]

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

                    centroid = np.mean(centers*resolutions, axis=0)

                    hyperstack, metadata = tif_reader_5D(path_data)
                    channels_seg = np.array([ch_A12, ch_F3, ch_Casp3])
                    hyperstack_seg = np.sum(hyperstack[:,:,channels_seg, :, :].astype("int32"), axis=2)

                    z_plot = np.rint(hyperstack_seg.shape[1]/2).astype("int64")
                    
                    if TIME=="96hr": 
                        binth = binths[f]
                    else:
                        binth=7.5
                    
                    print(binth)
                    ES = EmbryoSegmentation(
                            hyperstack_seg,
                            ksize=5,
                            ksigma=20,
                            binths=binth,
                            apply_biths_to_zrange_only=False,
                            checkerboard_size=10,
                            num_inter=100,
                            smoothing=20,
                            trange=None,
                            zrange=range(minz, maxz+1),
                            mp_threads=14,
                        )

                    ES(hyperstack_seg)
                    # ES.plot_segmentation(0, minz + 4)
                    # ES.plot_segmentation(0, z_plot)
                    # ES.plot_segmentation(0, maxz - 4)

        
                    import numpy as np
                    from skimage import measure

                    contour_points3D = []
                    for zid, z in enumerate(range(minz, maxz+1)):
                        contours = measure.find_contours(ES.LS[0][zid], 0.5)
                        contour = []
                        # Select the largest contour as the gastruloid contour
                        for cont in contours:
                            if len(cont)>len(contour):
                                contour = cont
                        for p in contour:
                            contour_points3D.append(np.array([z, p[1], p[0]]))
                    contour_points3D = np.array(contour_points3D)
                    
                    centers_Casp3 = []
                    centers_Casp3_F3 = []
                    centers_Casp3_A12 = []
                    for cell in CT_Casp3.jitcells:
                        centers_Casp3.append(cell.centers[0])
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
                        if mdiff > 0:
                            a12[idx] += mdiff
                        else: 
                            f3[idx] -= mdiff
        
                        if f3[idx] > a12[idx]:
                            centers_Casp3_F3.append(cell.centers[0])
                        else:
                            centers_Casp3_A12.append(cell.centers[0])

                    centers_Casp3 = np.array(centers_Casp3)
                    centers_Casp3_F3 = np.array(centers_Casp3_F3)
                    centers_Casp3_A12 = np.array(centers_Casp3_A12)

                    centers_A12 = []
                    for cell in CT_A12.jitcells:
                        centers_A12.append(cell.centers[0])
                    centers_A12 = np.array(centers_A12)

                    centers_F3 = []
                    for cell in CT_F3.jitcells:
                        centers_F3.append(cell.centers[0])
                    centers_F3 = np.array(centers_F3)

                    contour_points3D = contour_points3D*resolutions
                    
                    if len(centers_Casp3)!= 0:
                        centers_Casp3 = centers_Casp3*resolutions
                    if len(centers_Casp3_F3)!= 0:
                        centers_Casp3_F3 = centers_Casp3_F3*resolutions
                    if len(centers_Casp3_A12)!= 0:
                        centers_Casp3_A12 = centers_Casp3_A12*resolutions

                    centers_A12 = centers_A12*resolutions
                    centers_F3 = centers_F3*resolutions

                    if len(centers_Casp3)!=0:
                        dists_Casp3 = compute_dists_jit(np.array(centers_Casp3), np.array(contour_points3D), compute_distance_xyz_jit)
                        closests_Casp3_ids = np.argmin(dists_Casp3, axis=1)
                        closests_contour_points_Casp3 = np.array([contour_points3D[i] for i in closests_Casp3_ids])
                    else: 
                        closests_Casp3_ids = np.array([])
                        
                    if len(centers_Casp3_F3)!=0:
                        dists_Casp3_F3 = compute_dists_jit(np.array(centers_Casp3_F3), np.array(contour_points3D), compute_distance_xyz_jit)
                        closests_Casp3_ids_F3 = np.argmin(dists_Casp3_F3, axis=1)
                        closests_contour_points_Casp3_F3 = np.array([contour_points3D[i] for i in closests_Casp3_ids_F3])
                    else: 
                        closests_Casp3_ids_F3 = np.array([])
                        
                    if len(centers_Casp3_A12)!=0:
                        dists_Casp3_A12 = compute_dists_jit(np.array(centers_Casp3_A12), np.array(contour_points3D), compute_distance_xyz_jit)
                        closests_Casp3_ids_A12 = np.argmin(dists_Casp3_A12, axis=1)
                        closests_contour_points_Casp3_A12 = np.array([contour_points3D[i] for i in closests_Casp3_ids_A12])
                    else: 
                        closests_Casp3_ids_A12 = np.array([])
                        
                    dists_A12 = compute_dists_jit(np.array(centers_A12), np.array(contour_points3D), compute_distance_xyz_jit)
                    closests_A12_ids = np.argmin(dists_A12, axis=1)
                    closests_contour_points_A12 = np.array([contour_points3D[i] for i in closests_A12_ids])

                    dists_F3 = compute_dists_jit(np.array(centers_F3), np.array(contour_points3D), compute_distance_xyz_jit)
                    closests_F3_ids = np.argmin(dists_F3, axis=1)
                    closests_contour_points_F3 = np.array([contour_points3D[i] for i in closests_F3_ids])

                    dists_contour_Casp3_current = [dists_Casp3[i, closests_Casp3_ids[i]] for i in range(len(centers_Casp3))]
                    dists_contour_Casp3_current_F3 = [dists_Casp3_F3[i, closests_Casp3_ids_F3[i]] for i in range(len(centers_Casp3_F3))]
                    dists_contour_Casp3_current_A12 = [dists_Casp3_A12[i, closests_Casp3_ids_A12[i]] for i in range(len(centers_Casp3_A12))]

                    dists_contour_A12_current = [dists_A12[i, closests_A12_ids[i]] for i in range(len(centers_A12))]
                    dists_contour_F3_current = [dists_F3[i, closests_F3_ids[i]] for i in range(len(centers_F3))]
                    
                    dists_centroid_Casp3_current = [compute_distance_xyz_jit(center, centroid) for center in centers_Casp3]
                    dists_centroid_Casp3_current_F3 = [compute_distance_xyz_jit(center, centroid) for center in centers_Casp3_F3]
                    dists_centroid_Casp3_current_A12 = [compute_distance_xyz_jit(center, centroid) for center in centers_Casp3_A12]

                    dists_centroid_A12_current = [compute_distance_xyz_jit(center, centroid) for center in centers_A12]
                    dists_centroid_F3_current = [compute_distance_xyz_jit(center, centroid) for center in centers_F3]

                    file_path = path_save_results+embcode
                    np.save(file_path+"_dists_contour_Casp3", dists_contour_Casp3_current, allow_pickle=False)
                    np.save(file_path+"_dists_contour_Casp3_F3", dists_contour_Casp3_current_F3, allow_pickle=False)
                    np.save(file_path+"_dists_contour_Casp3_A12", dists_contour_Casp3_current_A12, allow_pickle=False)

                    np.save(file_path+"_dists_contour_A12", dists_contour_A12_current, allow_pickle=False)
                    np.save(file_path+"_dists_contour_F3", dists_contour_F3_current, allow_pickle=False)
                    
                    np.save(file_path+"_dists_centroid_Casp3", dists_centroid_Casp3_current, allow_pickle=False)
                    np.save(file_path+"_dists_centroid_Casp3_F3", dists_centroid_Casp3_current_F3, allow_pickle=False)
                    np.save(file_path+"_dists_centroid_Casp3_A12", dists_centroid_Casp3_current_A12, allow_pickle=False)

                    np.save(file_path+"_dists_centroid_A12", dists_centroid_A12_current, allow_pickle=False)
                    np.save(file_path+"_dists_centroid_F3", dists_centroid_F3_current, allow_pickle=False)
