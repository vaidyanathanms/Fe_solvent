#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 23:30:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --mem=2G
#SBATCH -J FeTFSI_0.4
#SBATCH -o out.%J
#SBATCH -e err.%J

module purge
module load gromacs/5.1.2


cd $SLURM_SUBMIT_DIR

echo "begin job.."
echo $PWD

mkdir -p outdir
mkdir -p trajfiles

#---------------------------------------------------------Generate initial run----------------------------------
finit_inp=./enermin.tpr
if ! test -f "$finit_inp"; then
	echo "begin generating enermin.tpr.."
	# generate enermin files
	srun gmx grompp -f minim.mdp -c main.pdb -p topol.top -o enermin.tpr
	
fi
wait

#---------------------------------------------------------Minimization -DFLEXIBLE----------------------------------
finit2_inp=./min_noflex.tpr
if ! test -f "$finit2_inp"; then
	
	echo "begin running enermin.tpr.."
	# run enermin.tpr
	srun gmx mdrun -s enermin.tpr -cpo state_min.cpt -cpi state_min.cpt -cpt 2 -g md_min.log -o traj_min.trr -e ener_min.edr -c confout_min.gro  
	wait

	echo "begin generating min_noflex.tpr.."
	# generate enermin files with no -DFLEXIBLE
	srun gmx grompp -f min_noflex.mdp -c confout_min.gro -p topol.top -o min_noflex.tpr
	wait

        cp md_min.log trajfiles/md_min.log
	cp traj_min.trr trajfiles/traj_min.trr
	cp ener_min.edr trajfiles/ener_min.edr
	cp confout_min.gro trajfiles/confout_min.gro
	
fi
wait

#----------------------------------------------------------Minimization no FLEXIBLE--------------------------------------
fmin_inp=./nvt.tpr 
if ! test -f "$fmin_inp"; then

	echo "begin running min_noflex.tpr.."
	# run enermin.tpr
	srun gmx mdrun -s min_noflex.tpr -cpo state_min2.cpt -cpi state_min2.cpt -cpt 2 -g md_min2.log -o traj_min2.trr -e ener_min2.edr -c confout_min2.gro  
	wait

	echo "begin generating nvt.tpr.."
	# generate nvt files
	srun gmx grompp -f nvt.mdp -c confout_min2.gro -p topol.top -o nvt.tpr 
	wait

        cp md_min.log trajfiles/md_min2.log
	cp traj_min.trr trajfiles/traj_min2.trr
	cp ener_min.edr trajfiles/ener_min2.edr
	cp confout_min.gro trajfiles/confout_min2.gro
fi
wait

#----------------------------------------------------------NVT Equilibration--------------------------------------------
fnvt_inp=./npt_berendsen.tpr
if ! test -f "$fnvt_inp"; then

	echo "begin running nvt.tpr.."
	# run nvt.tpr
	srun gmx mdrun -s nvt.tpr -cpo state_nvt.cpt -cpi state_nvt.cpt -cpt 5 -g md_nvt.log -o traj_nvt.trr -e ener_nvt.edr -c confout_nvt.gro  
	wait

	echo "begin generating npt_berendsen.tpr.."
	# generate npt_berendsen files
	srun gmx grompp -f npt_berendsen.mdp -c confout_nvt.gro -p topol.top -o npt_berendsen.tpr
	wait

	cp md_nvt.log trajfiles/md_nvt.log
	cp traj_nvt.trr trajfiles/traj_nvt.trr
	cp ener_nvt.edr trajfiles/ener_nvt.edr
	cp confout_nvt.gro trajfiles/confout_nvt.gro
fi
wait

#----------------------------------------------------------NPT Equilibration/Production--------------------------------------------
fnpt_inp=./npt_main.tpr
if ! test -f "$fnpt_inp"; then
         
 
        fnpt_cpt=./state_npt_berendsen.cpt 
        if test -f "$fnpt_cpt"; then
		echo "begin generating npt_berendsen.tpr from last npt state.."
		# generate npt_berendsen files
                srun gmx trjconv -f state_npt_berendsen.cpt -s npt_berendsen.tpr -o confout_last_npt.gro
                wait
		# Recreate npt_berendsen.tpr
		srun gmx grompp -f npt_berendsen.mdp -c confout_last_npt.gro -p topol.top -o npt_berendsen.tpr
		wait
	fi
	wait

	echo "begin running npt_berendsen.tpr.."
	# run npt_berendsen.tpr
	srun gmx mdrun -s npt_berendsen.tpr -cpo state_npt_berendsen.cpt -cpi state_npt_berendsen.cpt -cpt 5 -g md_npt_berendsen.log -o traj_npt_berendsen.trr -e ener_npt_berendsen.edr -c confout_npt_berendsen.gro  
	wait

	echo "begin generating npt_main.tpr.."
	# generate npt_berendsen files
	srun gmx grompp -f npt_main.mdp -c confout_npt_berendsen.gro -p topol.top -o npt_main.tpr
	wait

	cp md_npt_berendsen.log trajfiles/md_npt_berendsen.log
	cp traj_npt_berendsen.trr trajfiles/traj_npt_berendsen.trr
	cp ener_npt_berendsen.edr trajfiles/ener_npt_berendsen.edr
	cp confout_npt_berendsen.gro trajfiles/confout_npt_berendsen.gro
else

        echo "begin running npt_main.tpr.."
        # run npt_main.tpr
        srun gmx mdrun -s npt_main.tpr -cpo state_npt_main.cpt -cpi state_npt_main.cpt -cpt 5 -g md_npt_main.log -o traj_npt_main.trr -e ener_npt_main.edr -c confout_npt_main.gro  
        wait


        cp md_npt_main.log trajfiles/md_npt_main.log
        cp traj_npt_main.trr trajfiles/traj_npt_main.trr
        cp ener_npt_main.edr trajfiles/ener_npt_main.edr
        cp confout_npt_main.gro trajfiles/confout_npt_main.gro
fi
wait
#------------------------------------------------------------------------------------------------------------------------------------

echo "End of run.."
