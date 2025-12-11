### LOAD PACKAGE ###
from qlivecell import get_file_name, cellSegTrack, check_or_create_dir, get_file_names, correct_path, fill_channels

### LOAD STARDIST MODEL ###
from stardist.models import StarDist2D
model = StarDist2D.from_pretrained('2D_versatile_fluo')

EXPERIMENTS = ["2023_11_17_Casp3", "2024_03_Casp3"]
CONDS = ["WT", "KO"]

batch_args = {
    'name_format':"ch"+str(0)+"_{}",
    'extension':".tif",
} 

master_path_to_data = "/home/pablo/Desktop/PhD/projects/Data/gastruloids/joshi/competition/"
master_path_to_save = "/home/pablo/Desktop/papers/GastruloidCompetition_paper/cell_competition_gastruloids/image-analysis/results/segmentation/"

EXP = EXPERIMENTS[0]
if EXP=="2023_11_17_Casp3":
    TIMES = ["48hr", "72hr", "96hr"]
elif EXP=="2024_03_Casp3":        
    TIMES = ["48hr", "60hr", "72hr", "96hr"]
    
TIME = TIMES[0]

COND = CONDS[0]
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

f = 0
file = files[f]
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
CT_F3.plot_tracking()
    
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
CT_A12.plot_tracking()
