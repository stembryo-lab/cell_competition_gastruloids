### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path
import numpy as np
import matplotlib.pyplot as plt

### LOAD STARDIST MODEL ###
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained('2D_versatile_fluo')

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
CONDS = ["WT", "KO"]

### DEFINE ARGUMENTS ###
segmentation_args={
    'method': 'stardist2D', 
    'model': model, 
    'blur': [2,1], 
    'min_outline_length':100,
}

concatenation3D_args = {
    'distance_th_z': 3.0, # microns
    'relative_overlap':False, 
    'use_full_matrix_to_compute_overlap':True, 
    'z_neighborhood':2, 
    'overlap_gradient_th':0.3, 
    'min_cell_planes': 2,
}

error_correction_args = {
    'backup_steps': 10,
    'line_builder_mode': 'points',
}
 
master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"
               
for EXP in EXPERIMENTS:
    if EXP=="2023_11_17_Casp3":
        TIMES = ["48hr", "72hr", "96hr"]
    elif EXP=="2024_03_Casp3":        
        TIMES = ["48hr", "60hr", "72hr", "96hr"]
        
    for TIME in TIMES:
        for COND in CONDS:
            path_data_dir=correct_path(master_path_to_data)+"{}/stacks/{}/{}/".format(EXP, TIME, COND)
            path_save_dir=correct_path(master_path_to_save)+"{}/segmentation_results/{}/{}/".format(EXP, TIME, COND)
            check_or_create_dir(path_save_dir)

            files = get_file_names(path_data_dir)
            
            if EXP=="2023_11_17_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
                if "96hr" in path_data_dir:
                    channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
            elif EXP=="2024_03_Casp3":
                channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]

            for f, file in enumerate(files):
                path_data = path_data_dir+file
                file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
                path_save = path_save_dir+embcode
                check_or_create_dir(path_save)


                
                ch = channel_names.index("F3")
                
                batch_args = {
                    'name_format':"ch"+str(ch)+"_{}",
                    'extension':".tif",
                } 
                plot_args['channels'] = [ch]
                
                chans = [ch]
                for _ch in range(len(channel_names)):
                    if _ch not in chans:
                        chans.append(_ch)
                        
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
                # CT_F3.plot(plot_args=plot_args)
                    
                ch = channel_names.index("A12")
                batch_args = {
                    'name_format':"ch"+str(ch)+"_{}",
                    'extension':".tif",
                }
                plot_args = {
                    'plot_layout': (1,1),
                    'plot_overlap': 1,
                    'masks_cmap': 'tab10',
                    'plot_stack_dims': (512, 512), 
                    'plot_centers':[False, False], # [Plot center as a dot, plot label on 3D center]
                    'channels':[ch],
                    'min_outline_length':75,
                }

                chans = [ch]
                for _ch in range(len(channel_names)):
                    if _ch not in chans:
                        chans.append(_ch)

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
                # CT_A12.plot(plot_args=plot_args)