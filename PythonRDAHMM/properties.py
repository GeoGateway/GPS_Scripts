#!/usr/local/bin/python
#==========================================================================
# Definitions of all global variables such as paths and commands used by
# data processing scripts. Imported and invoked internally. 
#
#===========================================================================

# def properties(key):
    # V={}
    # V['cron_path']="/home/yuma/RDAHMM/CRON_Download/"  
    # V['download_path']="/home/yuma/RDAHMM/Download/"  
    # V['script_path']="/home/yuma/PythonRDAHMM/"
    # V['data_path']="/home/yuma/RDAHMM/Data/"
    # # temp_path is the temporary working directory for ingesting raw data
    # V['temp_path']="/home/yuma/RDAHMM/TEMP/"
    # V['model_path']="/home/yuma/RDAHMM/Model/"
    # V['eval_path']="/var/www/html/daily_rdahmmexec/daily/"
    # V['train_epoch']="2013-12-31"
    # V['rdahmm_bin']="/home/yuma/RDAHMM/rdahmm3/bin/rdahmm"
    # V['rdahmm_model_parm']="-data <inputFile> -T <dataCount> -D <dimensionCount> -N 5 -output_type gauss -anneal -annealfactor 1.1 -betamin 0.1 -regularize -omega 0 0 1 1.0e-6 -ntries 10 -seed 1234"
    # V['rdahmm_eval_parm']="-data <proBaseName>.all.input -T <dataCount> -D <dimensionCount> -N 5 -output_type gauss -A <modelBaseName>.A -B <modelBaseName>.B -pi <modelBaseName>.pi -minvalfile <modelBaseName>.minval -maxvalfile <modelBaseName>.maxval -rangefile <modelBaseName>.range -eval"
    # V['dygraphsJs']="/home/yuma/PythonRDAHMM/dygraphsJsCreator.perl"
    # return V[key]

import os
def get_parent_dir(directory):
# http://stackoverflow.com/questions/9856683/using-pythons-os-path-how-do-i-go-up-one-directory
    return os.path.dirname(directory)

Paths={}
# WORK_DIR = get_parent_dir(os.getcwd())
# WORK_DIR = os.getcwd()
WORK_DIR = '~/GPS_Scripts'
Paths['cron_path']=os.path.join(WORK_DIR, "RDAHMM","CRON_Download/")  
Paths['download_path']=os.path.join(WORK_DIR,"RDAHMM","Download/")  
Paths['script_path']=os.path.join(WORK_DIR,"PythonRDAHMM/") 
Paths['data_path']=os.path.join(WORK_DIR,"RDAHMM","Data/")
# temp_path is the temporary working directory for ingesting raw data
Paths['temp_path']=os.path.join(WORK_DIR,"RDAHMM","TEMP/")
Paths['model_path']=os.path.join(WORK_DIR,"RDAHMM","Model/")
Paths['eval_path']=os.path.join(WORK_DIR,"daily/")
Paths['train_epoch']="2013-12-31"
Paths['rdahmm_bin']=os.path.join(WORK_DIR,"RDAHMM", "rdahmm3","bin", "rdahmm")
Paths['rdahmm_model_parm']="-data <inputFile> -T <dataCount> -D <dimensionCount> -N 5 -output_type gauss -anneal -annealfactor 1.1 -betamin 0.1 -regularize -omega 0 0 1 1.0e-6 -ntries 10 -seed 1234"
Paths['rdahmm_eval_parm']="-data <proBaseName>.all.input -T <dataCount> -D <dimensionCount> -N 5 -output_type gauss -A <modelBaseName>.A -B <modelBaseName>.B -pi <modelBaseName>.pi -minvalfile <modelBaseName>.minval -maxvalfile <modelBaseName>.maxval -rangefile <modelBaseName>.range -eval"
Paths['dygraphsJs']=os.path.join(WORK_DIR,"PythonRDAHMM","dygraphsJsCreator.pearl")

def properties(key):
    return Paths[key]