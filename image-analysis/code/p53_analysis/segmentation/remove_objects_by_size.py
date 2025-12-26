### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, fill_channels

# ### PATH TO YOU DATA FOLDER AND TO YOUR SAVING FOLDER ###
CONDS = ["WT", "KO"]
repeats = ["n2", "n3", "n4", "n3_Ab"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/segmentation_results/"
path_to_save_results = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/p53/figures/"
check_or_create_dir(path_to_save_results)

size_th_low = 31.9
size_th_high = 250.0
for COND in CONDS:
    for REP in repeats:
        
        path_data_dir="/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/p53_analysis/input/{}/{}/".format(COND,REP)
        path_save_dir = "{}{}/{}/".format(master_path_save, COND, REP)
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
            chans = fill_channels(channel=ch, channel_names=channel_names)
            
            CT_F3 = cellSegTrack(
                path_data,
                path_save,
                batch_args=batch_args,
                channels=chans
            )

            CT_F3.load()
            
            labs_to_rem = []
            for cell in CT_F3.jitcells:
                zc = int(cell.centers[0][0])
                zcid = cell.zs[0].index(zc)

                mask = cell.masks[0][zcid]
                area = len(mask) * CT_F3.metadata["XYresolution"]**2
                if area < size_th_low:
                    labs_to_rem.append(cell.label)
                elif area > size_th_high:
                    labs_to_rem.append(cell.label)
                    
            for lab in labs_to_rem:
                CT_F3._del_cell(lab)  
                 
            CT_F3.update_labels()
            
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
            for cell in CT_A12.jitcells:
                zc = int(cell.centers[0][0])
                zcid = cell.zs[0].index(zc)

                mask = cell.masks[0][zcid]
                area = len(mask) * CT_A12.metadata["XYresolution"]**2
                if area < size_th_low:
                    labs_to_rem.append(cell.label)
                elif area > size_th_high:
                    labs_to_rem.append(cell.label)
            for lab in labs_to_rem:
                CT_A12._del_cell(lab)  
                 
            CT_A12.update_labels()
