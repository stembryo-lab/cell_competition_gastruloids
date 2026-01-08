"""
Main pipeline for 3D cell segmentation of gastruloid imaging data.

This script performs:
- 3D cell segmentation using a StarDist2D model via the qlivecell framework
- Concatenation of 2D segmentations across z-planes to obtain 3D cell objects
- Independent segmentation of mCherry and emiRFP channels

This script is intended solely as a preprocessing step to generate 3D segmentation
results. Quantitative analysis (e.g. cell counts, density, neighborhood composition,
or radial distribution) is not performed in this file.
"""

### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels

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

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"

for EXP in EXPERIMENTS:
    if EXP=="2023_11_17_Casp3":
        TIMES = ["48hr", "72hr", "96hr"]
    elif EXP=="2024_03_Casp3":        
        TIMES = ["48hr", "60hr", "72hr", "96hr"]
        
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
            check_or_create_dir(path_save_dir)
            files = get_file_names(path_data_dir)
            
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
                    segmentation_args=segmentation_args,
                    concatenation3D_args=concatenation3D_args,
                    batch_args=batch_args,
                    channels=chans
                )

                CT_F3.run()

                ch = channel_names.index("A12")
                batch_args['name_format'] = "ch"+str(ch)+"_{}"                

                chans = fill_channels(channel=ch, channel_names=channel_names)

                CT_A12 = cellSegTrack(
                    path_data,
                    path_save,
                    segmentation_args=segmentation_args,
                    concatenation3D_args=concatenation3D_args,
                    batch_args=batch_args,
                    channels=chans
                )

                CT_A12.run()
