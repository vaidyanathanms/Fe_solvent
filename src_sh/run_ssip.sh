#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 05:00:00
#SBATCH -N 1                    # 1 nodes
#SBATCH --ntasks-per-node=1     # 1 tasks/node
#SBATCH -c 1                    # 1 cores/task
#SBATCH --mem=0
#SBATCH -J ssip_Fe3_m_5.0
#SBATCH -o ana_mol_5.0/out.%J
#SBATCH -e ana_mol_5.0/err.%J

# Load modules
module reset
module load intel mkl
module load python/3.15.04

# Export gcc path and num_threads
export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
cd $SLURM_SUBMIT_DIR
echo "begin job to analyze SSIP  @start time: ${date}"
echo $PWD

outdir='ana_mol_5.0'
mkdir -p ${outdir}
wait

python ssip.py -s npt_main.tpr -f traj_npt_main_nojump.trr
wait

mv pairs_summary.csv ${outdir}
echo "End of run.."
