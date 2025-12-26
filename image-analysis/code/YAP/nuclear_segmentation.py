from qlivecell import get_file_name, cellSegTrack, get_file_names, check_or_create_dir, fill_channels
from cellpose import models
model = models.CellposeModel(gpu=True, model_type='cyto3')
channel_names = ["A12", "GPI-GFP", "YAP", "DAPI"]

size_th = 30.0 #µm

ch = channel_names.index("DAPI")

segmentation_args={
'method': 'cellpose2D', 
'model': model, 
'blur': None, 
'channels': [ch+1, 0],
'diameter': 200,
}

concatenation3D_args = {
    'do_3Dconcatenation':False
}

error_correction_args = {
    'backup_steps': 10,
    'line_builder_mode': 'points',
}

chans = fill_channels(ch, channel_names=channel_names)

batch_args = {
    'name_format':"ch"+str(ch)+"_{}",
    'extension':".tif",
}

plot_args = {
    'plot_layout': (1,1),
    'plot_overlap': 1,
    'masks_cmap': 'tab10',
    'plot_stack_dims': (256, 256), 
    'plot_centers':[False, False], # [Plot center as a dot, plot label on 3D center]
    'channels':[ch],
    'min_outline_length':75,
}

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/YAP/segmentation_results/"
check_or_create_dir(master_path_save)

CONDITIONS = ["WT", "KO8", "KO25", "KO8_ABonly", "KO25_ABonly"]
for COND in CONDITIONS:
    path_data_dir = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/YAP/2025_02_02_AiryscMultipl_FastMediumQuality_Files/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    
    ### GET FULL FILE NAME AND FILE CODE ###
    files = get_file_names(path_data_dir)
    
    for file in files:
        if not ".tif" in file: continue
        
        file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
        
        path_data = path_data_dir+file
        path_save = path_save_dir+embcode
        check_or_create_dir(path_save)
            
        CT = cellSegTrack(
            path_data,
            path_save,
            segmentation_args=segmentation_args,
            concatenation3D_args=concatenation3D_args,
            error_correction_args=error_correction_args,
            plot_args=plot_args,
            batch_args=batch_args,
            channels=chans
        )
        CT.run()
        
        labs_to_rem = []
        for cell in CT.jitcells:
            zc = int(cell.centers[0][0])
            zcid = cell.zs[0].index(zc)

            mask = cell.masks[0][zcid]
            area = len(mask) * CT.metadata["XYresolution"]**2
            if area < size_th:
                labs_to_rem.append(cell.label)
            
        for lab in labs_to_rem:
            CT._del_cell(lab)  

        CT.update_labels()
        # CT.plot_tracking(plot_args=plot_args)
