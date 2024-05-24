from pathlib import Path as path_

files = []
refs = []

ROOT_DIR = path_(__file__).parent
init_directory = ROOT_DIR

# =================================== bools ================================== #
# save_mean_bool = True
# save_std_time_bool = True
# save_std_freq_bool = True
save_bools = {
    "mean" : True,
    "std_time" : True,
    "std_freq" : True,
    "correction_param" : False,
    "time_traces" : False,
    "noise_matrix" : False,
}

myinput=[]
myreferencedata = []
ncm = []
ncm_inverse = []
reference_number = []
mydatacorrection = []
delay_correction = []
dilatation_correction = []
leftover_correction = []
myglobalparameters = []
fopt = []
fopt_init = []
mode = []
preview = []

dialog_initialzed = False
modesuper = None


