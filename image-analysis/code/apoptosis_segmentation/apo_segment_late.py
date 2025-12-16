### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels
import numpy as np

### LOAD STARDIST MODEL ###
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained('2D_versatile_fluo')

apo_stage = "late"
CONDS = ["WT", "KO"]

### DEFINE ARGUMENTS ###
segmentation_args={
    'method': 'stardist2D', 
    'model': model, 
    'blur': [1,1], 
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
    'line_builder_mode': 'points',
}
             
plot_args = {
    'plot_layout': (1,1),
    'plot_overlap': 1,
    'masks_cmap': 'tab10',
    'plot_stack_dims': (512, 512), 
    'plot_centers':[False, False], # [Plot center as a dot, plot label on 3D center]
    'channels':[0],
    'min_outline_length':75,
}
                   
batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

for experiment in ["2023_11_17_Casp3", "2024_03_Casp3"]:

    if experiment=="2023_11_17_Casp3":
        channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
        if "96hr" in path_data_dir:
            channel_names = ["A12", "F3", "Casp3", "BF", "DAPI"]
        TIMES = ["48hr", "72hr", "96hr"]
    else: 
        channel_names = ["F3", "A12", "DAPI", "Casp3", "BF"]
        TIMES = ["48hr", "60hr", "72hr", "96hr"]

    for TIME in TIMES:
        for COND in CONDS:
            path_data_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/stacks/{}/{}/'.format(experiment, TIME, COND)
            path_save_dir='/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/{}/ctobjects/{}/{}/'.format(experiment, TIME, COND)

            check_or_create_dir(path_save_dir)
            files = get_file_names(path_data_dir)
        
            for f, file in enumerate(files):
                
                path_data = path_data_dir+file
                file, embcode = get_file_name(path_data_dir, file, allow_file_fragment=False, return_files=False, return_name=True)
                path_save = path_save_dir+embcode
                check_or_create_dir(path_save)

                ch_F3 = channel_names.index("F3")
                ch_A12 = channel_names.index("A12")
                ch_Casp3 = channel_names.index("Casp3")
                ch_DAPI = channel_names.index("DAPI")
            
                chans = fill_channels(channel=ch_Casp3, channel_names=channel_names)
                  
                batch_args['name_format'] = "ch"+str(ch_Casp3)+"_{}_"+apo_stage                
                plot_args['channels'] = [ch_A12, ch_F3, ch_Casp3]
                
                CT_Casp3 = cellSegTrack(
                    path_data,
                    path_save,
                    segmentation_args=segmentation_args,
                    error_correction_args=error_correction_args,
                    plot_args=plot_args,
                    batch_args=batch_args,
                    channels=chans
                )

                CT_Casp3.run()
                        
                intens = []
                for cell in CT_Casp3.jitcells:
                    zc = int(cell.centers[0][0])
                    zcid = cell.zs[0].index(zc)
                    mask = cell.masks[0][zcid]
                    ints = CT_Casp3.hyperstack[0,zc,ch_Casp3,:,:][mask[:,1], mask[:,0]]
                    intens = [*intens, *ints]
                
                int_th = np.percentile(intens, 90)

                labs_to_rem = []
                for cell in CT_Casp3.jitcells:
                    zc = int(cell.centers[0][0])
                    zcid = cell.zs[0].index(zc)
                    mask = cell.masks[0][zcid]
                    area = len(mask)* CT_Casp3.CT_info.xyresolution**2
                    ints = np.mean(CT_Casp3.hyperstack[0,zc,ch_Casp3,:,:][mask[:,1], mask[:,0]])
                    if area > 31.75:
                        labs_to_rem.append(cell.label)
                    if ints < int_th and cell.label not in labs_to_rem:
                        labs_to_rem.append(cell.label)
                        
                for lab in labs_to_rem:
                    CT_Casp3._del_cell(lab)   
                