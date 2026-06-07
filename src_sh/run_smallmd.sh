#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 12:00:00
#SBATCH -N 4                    # 4 nodes
#SBATCH --ntasks-per-node=4     # 4 tasks/node
#SBATCH -c 8                    # 8 cores/task
#SBATCH --mem=0
#SBATCH -J FeTFSI_2.5
#SBATCH -o /home/vm5/outdir/out.%J
#SBATCH -e /home/vm5/outdir/err.%J

# Load modules
module reset
module load intel mkl fftw hwloc cmake

# Export gcc path and num_threads
export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
cd /lustre/or-scratch/cades-birthright/vm5/fetfsi/fetfsi_3/trial_4/mol_2.5
echo "begin job @start time: ${date}"
echo $PWD

# Run command
gmx="${HOME}/gromacs-2024.5/install/bin/gmx_mpi "

mkdir -p outdir

echo "begin generating npt_smallrun.tpr.."
# generate npt_berendsen files
srun ${gmx} grompp -f npt_smallrun.mdp -c confout_npt_main.gro -p topol.top -o npt_smallrun.tpr
wait

echo "begin running npt_smallrun.tpr.."
# run npt_main.tpr
srun ${gmx} mdrun -s npt_smallrun.tpr -cpo state_npt_smallrun.cpt -cpi state_npt_smallrun.cpt -cpt 5 -g md_npt_smallrun.log -o traj_npt_smallrun.trr -e ener_npt_smallrun.edr -c confout_npt_smallrun.gro -notunepme 
wait
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."
