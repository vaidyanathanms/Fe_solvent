#!/bin/bash

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=2G
#SBATCH -J rdfana_0.1
#SBATCH -o rdfout.%J
#SBATCH -e rdferr.%J

module purge
module load gromacs/5.1.2


cd $SLURM_SUBMIT_DIR

echo "begin job.."
echo $PWD

# Make output directory
outdir='rdfresults_0.1'
mkdir -p ${outdir}

# Create RDF files
#echo "Run rdf-Wat_Wat"
#srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Wat_Wat.xvg -cn nrdf_Wat_Wat.xvg   -b 80000 -e 90000 -rmpbc yes -sf rdfWat_Wat.txt
#wait

#echo "Run rdf-TFSI-TFSI"
#srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_TFSI_TFSI.xvg -cn nrdf_TFSI_TFSI.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfTFSI_TFSI.txt
#wait

#echo "Run rdf-Wat_TFSI"
#srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Wat_TFSI.xvg -cn nrdf_Wat_TFSI.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfWat_TFSI.txt
#wait

echo "Run rdf-Fe-TFSI"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_TFSI.xvg -cn nrdf_Fe_TFSI.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_TFSI.txt
wait

#echo "Run rdf-Fe_Wat"
#srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_Wat.xvg -cn nrdf_Fe_Wat.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_Wat.txt
#wait

#echo "Run rdf-Fe-Fe"
#srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_Fe.xvg -cn nrdf_Fe_Fe.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_Fe.txt
#wait

echo "Run rdf-Fe-OTFSI"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_OTFSI.xvg -cn nrdf_Fe_OTFSI.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_OTFSI.txt
wait

echo "Run rdf-Fe-NTFSI"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_NTFSI.xvg -cn nrdf_Fe_NTFSI.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_NTFSI.txt
wait

echo "Run rdf-Fe-OWat"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_OWat.xvg -cn nrdf_Fe_OWat.xvg -b 80000 -e 90000 -rmpbc yes -sf rdfFe_owat.txt
wait


echo "All RDF calculations completed"
echo "move files to ${outdir}"
mv rdf_*xvg ${outdir} 
cp rdf*txt ${outdir}
