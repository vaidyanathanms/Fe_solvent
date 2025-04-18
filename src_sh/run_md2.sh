#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 23:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=2G
#SBATCH -J FeTFSI_0.1
#SBATCH -o out.%J
#SBATCH -e err.%J

module load openmpi
module load python
module load gromacs


cd $SLURM_SUBMIT_DIR

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

#---------------------------------------------------------Generate initial run----------------------------------
fnpt_inp=./npt_main2.tpr
if ! test -f "$fnpt_inp"; then
	echo "begin generating npt_main2.tpr.."
	# generate enermin files
	srun gmx_mpi grompp -f npt_main2.mdp -c confout_npt_main.gro -p topol.top -o npt_main2.tpr
	
fi
wait

echo "begin running npt_main2.tpr.."
# run npt_main.tpr
srun gmx_mpi mdrun -s npt_main2.tpr -cpo state_npt_main2.cpt -cpi state_npt_main2.cpt -cpt 5 -g md_npt_main2.log -o traj_npt_main2.trr -e ener_npt_main2.edr -c confout_npt_main2.gro   
wait

cp md_npt_main2.log trajfiles/md_npt_main2.log
cp traj_npt_main2.trr trajfiles/traj_npt_main2.trr
cp ener_npt_main2.edr trajfiles/ener_npt_main2.edr
cp confout_npt_main2.gro trajfiles/confout_npt_main2.gro
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."
