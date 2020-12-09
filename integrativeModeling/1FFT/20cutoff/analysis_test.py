import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('/home/cns-mccafferty/PMI_analysis/pyext/src/')
from analysis_trajectories import *

current_dir = '/home/cns-mccafferty/Ecoli_CoEv_modeling/CoEv_IMP/1fft/20cutoff/'

# How are the trajectories dir names
dir_head = 'run'
out_dirs = glob.glob(current_dir+'/'+dir_head+'*/output/')

AT = AnalysisTrajectories(out_dirs,
                          dir_name = dir_head,
                          analysis_dir = current_dir + 'analysis',
                          nproc=16)
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
AT.set_analyze_Distance_restraint()

# Read stat files
AT.read_stat_files()
AT.write_models_info()
AT.hdbscan_clustering(['EV_sum', 'CR_sum', 'DR_sum'], min_cluster_size=500)
