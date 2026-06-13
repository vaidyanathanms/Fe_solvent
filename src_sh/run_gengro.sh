#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 01:00:00
#SBATCH -N 1                    # 1 nodes
#SBATCH --ntasks-per-node=1     # 1 tasks/node
#SBATCH -c 8                    # 8 cores/task
#SBATCH --mem=0
#SBATCH -J gengro_Fe3_m_0.1
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

# Load modules
module reset
module load intel mkl fftw hwloc cmake

# Export gcc path and num_threads
export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
cd $SLURM_SUBMIT_DIR
echo "begin job to generate gro files for analyses @start time: ${date}"
echo $PWD

# Run command
gmx="${HOME}/gromacs-2024.5/install/bin/gmx_mpi "

# Compute diffusivity

# write gro file at the msd start time - have to use the FULL system and not parts of it
echo "Generating confout file at 60 ps"
echo "0" | srun ${gmx} trjconv -f traj_npt_main.trr -s npt_main.tpr -o confout_all_60ps.gro -pbc whole -b 60000 -e 60000
wait

# generate tpr file at the msd start time
echo "Generating tpr file at 60 ps"
srun ${gmx} grompp -f npt_main.mdp -c confout_all_60ps.gro -p topol.top -o npt_main_60ps.tpr
wait

# write gro file at the t=0 have to use the FULL system and not parts of it
echo "Generating confout file at t=0 in NPT cycle"
echo "0" | srun ${gmx} trjconv -f traj_npt_main.trr -s npt_main.tpr -o confout_all_ps.gro -pbc whole -b 0 -e 0
wait

# convert npt_main.trr with no jump for the region of interest - here use the first frame and not the frame at 60 ps
echo "Convert to pbc nojump for the regions of interest"
echo "0" | srun ${gmx} trjconv -f traj_npt_main.trr -s npt_main.tpr -o traj_npt_main_nojump.gro -pbc nojump -dt 10 -b 50000 -e 99900 
wait


echo "End of run.."
