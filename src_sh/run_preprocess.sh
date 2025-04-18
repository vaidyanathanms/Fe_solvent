#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --mem=2G
#SBATCH -J FeTFSI_0.2
#SBATCH -o out.%J
#SBATCH -e err.%J

module purge
module load gromacs/5.1.2


cd $SLURM_SUBMIT_DIR

echo "begin job.."
echo $PWD

mkdir -p initdir
mkdir -p outdir

# editconf box
srun gmx editconf -f FeTFSI.pdb -bt cubic -box 6.40 6.40 6.40 -o main.pdb
wait

# make tpr file
srun gmx grompp -f minim.mdp -p topol.top -c main.pdb -o enermin.tpr
wait

cp *.pdb initdir/
