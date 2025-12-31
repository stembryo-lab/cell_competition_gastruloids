### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels

### LOAD STARDIST MODEL ###
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained('2D_versatile_fluo')

### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
CONDS = ["auxin_48-72_48", "auxin_48-72_72" , "auxin_48-72_96", "auxin_72-96_72", "auxin_72-96_96", "noauxin_72", "noauxin_96", "secondaryonly"]

### DEFINE ARGUMENTS ###
segmentation_args={
    'method': 'stardist2D', 
    'model': model, 
    'blur': [2,1], 
    'min_outline_length':100,
}

concatenation3D_args = {
    'do_3Dconcatenation': False
}

error_correction_args = {
    'backup_steps': 10,
    'line_builder_mode': 'points',
}

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 
plot_args = {
    'plot_layout': (1,1),
    'plot_overlap': 1,
    'masks_cmap': 'tab10',
    'plot_stack_dims': (256, 256), 
    'plot_centers':[False, False], # [Plot center as a dot, plot label on 3D center]
    'channels':[0],
    'min_outline_length':75,
}

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53_degron/segmentation_results/"
check_or_create_dir(master_path_save)

for COND in CONDS:        
    path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/2025_09_09_OsTIRMosaic_p53Timecourse/{}/".format(COND)
    path_save_dir = "{}{}/".format(master_path_save, COND)
    check_or_create_dir(path_save_dir)
    files = get_file_names(path_data_dir)
    
    channel_names = ["A12", "p53", "F3", "DAPI"]

    for f, file in enumerate(files):
        
        path_data = path_data_dir+file
        file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
        path_save = path_save_dir+embcode
        check_or_create_dir(path_save)

        ch = channel_names.index("F3")
        batch_args['name_format'] = "ch"+str(ch)+"_{}"    
        plot_args['channels'] = [ch]  
        chans = fill_channels(channel=ch, channel_names=channel_names)
                
        CT_F3 = cellSegTrack(
            path_data,
            path_save,
            segmentation_args=segmentation_args,
            concatenation3D_args=concatenation3D_args,
            error_correction_args=error_correction_args,
            plot_args=plot_args,
            batch_args=batch_args,
            channels=chans
        )

        CT_F3.run()
                    
        ch = channel_names.index("A12")
        batch_args['name_format'] = "ch"+str(ch)+"_{}"  
        plot_args['channels'] = [ch]    
        chans = fill_channels(channel=ch, channel_names=channel_names)

        CT_A12 = cellSegTrack(
            path_data,
            path_save,
            segmentation_args=segmentation_args,
            concatenation3D_args=concatenation3D_args,
            error_correction_args=error_correction_args,
            plot_args=plot_args,
            batch_args=batch_args,
            channels=chans
        )

        CT_A12.run()