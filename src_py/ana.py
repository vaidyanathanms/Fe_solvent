import numpy
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call


#---------mypython functions------------------------------

from my_python_functions import cpy_main_files
from my_python_functions import compile_anafiles_fe
from my_python_functions import find_trajfiles
from my_python_functions import edit_generate_anainp_files
from my_python_functions import run_analysis

#---------input details----------------------------------------
valency_fe   = '3' # 2 or 3 as string
trial_num    = '5' # trial number in the directory as string
analyze_only = 'all' #latest, all, filename, filelist
mol_fetfsi   = ['0.1','0.2','0.4','1.0']#,'2.5','5.0']# mol fraction of fetfsi
nframes      = 2000 # total frames to be analyzed
skipfr       = 100 # skip frames
start_time   = 52000 # for gromacs
freqfr       = 1 # freq of anaylsis
runana       = 1 # 1 - run analysis

#---------job details------------------------------------------
tottime   = 18 # in hours
nnodes    = 1 # number of nodes
ncores    = 8 # number of cores
hpc_sys   = 'cades'  # Opt: kestrel, cades

#--------file_lists--------------------------------------------
ana_files = ['analyze_fe.f90','gmx_params.f90','anainp_var.txt']
job_files = ['jobana_var.sh']
traj_pref = 'traj_npt_*nojump.gro'

#---------directory info---------------------------------------
maindir = os.getcwd() #src_py dir
if hpc_sys == 'kestrel':
    home_path = '/home/vaidyams'
    scr_path  = '/scratch/vaidyams'
elif hpc_sys == 'cades':
    home_path = '/home/vm5'
    scr_path  = '/lustre/or-scratch/cades-birthright/vm5'
else:
    raise RuntimeError('Unknown HPC system ' + hpc_sys)

src_f90 = home_path + '/all_codes/files_fesystems/src_f90' #src_f90 dir
src_sh  = home_path + '/all_codes/files_fesystems/src_sh' #src_f90 dir

scratchdir = scr_path  + '/fetfsi' #output headdir
scr_head   = 'fetfsi_' + str(valency_fe) # head dir scratch'

runjob_file = src_sh + '/runana_var.sh'
#--------lammps executable-------------------------------------
f90_comp   = 'ifx' # FORTRAN compiler

if not os.path.isdir(scratchdir):
    raise RuntimeError(scratchdir + " not found!")

#---------main analysis---------------------------------------
for molval in range(len(mol_fetfsi)):
    
    print(f"Analyzing Fetfsi_{valency_fe} for {mol_fetfsi[molval]} m")
    workdir1 = scratchdir + '/' + scr_head + '/trial_' + trial_num
    
    if not os.path.isdir(workdir1):
        print("ERROR: " + workdir1 + " not found!"); continue

    workdir_main = workdir1 + '/mol_' + str(mol_fetfsi[molval])
	
    if not os.path.isdir(workdir_main):
        print("ERROR: " + workdir_main + " not found!"); continue
        

    os.chdir(workdir_main)
    destdir = os.getcwd()

    #---Make a results director
    if not os.path.isdir('allresults_mol_'+str(mol_fetfsi[molval])):
        os.mkdir('allresults_mol_'+str(mol_fetfsi[molval]))

    #---Copying files------
    print( "Current Dir ", destdir)
    print( "Copying Files")
                
    for fyllist in range(len(ana_files)):
        cpy_main_files(src_f90,destdir,ana_files[fyllist])

    print( "Compiling analysis codes ...")
    compile_anafiles_fe(f90_comp)            

    #----Retrieve trajectory files
    print(" Finding trajectory files...")
    if analyze_only == 'filelist':
        if not os.path.exists(analist[fr_an]):
            print("ERROR: " + analist[fr_an] + " not found!")
            continue
        traj_arr = find_trajfiles(analyze_only,traj_pref,\
                                          analist[fr_an])
    else:
        traj_arr = find_trajfiles(analyze_only,traj_pref,\
                                          'none')
    if traj_arr == []:
        print("ERROR: No trajectory files found"); continue

    #----Iterate through trajectory files and submit-----            
    for fyllist in range(len(traj_arr)):
        print("Analyzing ",fyllist, traj_arr[fyllist])
        if 'smallrun' in traj_arr[fyllist]:
            dataname = 'confout_npt_main.gro'
            jobstr = 'ana_small_Fe' +str(valency_fe)+'_'+\
                str(mol_fetfsi[molval])
            jobana = 'runana_short.sh'
            start_time = 0
        elif 'main' in  traj_arr[fyllist]:
            dataname = 'confout_all_60ps.gro'
            jobstr = 'ana_main_Fe' +str(valency_fe)+'_'+\
                str(mol_fetfsi[molval])
            jobana = 'runana_main.sh'
            start_time = 52000
        else:
            print('Unknown dataname')
            continue
        anainp = edit_generate_anainp_files(dataname,traj_arr[fyllist],\
                                            0,nframes,skipfr,freqfr,\
                                            fyllist+1,start_time)
        print(anainp)
        run_analysis(anainp, jobstr,0,fyllist+1,\
                     runjob_file,jobana,traj_arr[fyllist],\
                     tottime,nnodes,ncores,runana)
        
