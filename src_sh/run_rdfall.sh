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
inpdir='rdfinps'
outdir='rdfall_0.1'
mkdir -p ${outdir}

# Create RDF files
echo "Run rdf-Fe-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_Fe_all.xvg -cn nrdf_Fe_all.xvg   -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refFe.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-Hwat-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_HWat_all.xvg -cn nrdf_HWat_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refHWat.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-Owat-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_OWat_all.xvg -cn nrdf_OWat_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refOWat.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-NTFSI-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_NTFSI_all.xvg -cn nrdf_NTFSI_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refNTFSI.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-OTFSI-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_OTFSI_all.xvg -cn nrdf_OTFSI_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refOTFSI.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-FTFSI-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_FTFSI_all.xvg -cn nrdf_FTFSI_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refFTFSI.txt -sel -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-CTFSI-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_CTFSI_all.xvg -cn nrdf_CTFSI_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refCTFSI.txt -sel  -sf ${inpdir}/rdfsel_all.txt
wait

echo "Run rdf-STFSI-all"
srun gmx rdf -f traj_npt_main.trr -s npt_main.tpr -o rdf_STFSI_all.xvg -cn nrdf_STFSI_all.xvg -b 80000 -e 90000 -rmpbc yes -ref -sf ${inpdir}/refSTFSI.txt -sel  -sf ${inpdir}/rdfsel_all.txt
wait


echo "All RDF calculations completed"
echo "move files to ${outdir}"
mv rdf_*xvg ${outdir} 
mv nrdf_*xvg ${outdir}
cp rdf*txt ${outdir}
rm *_Sofk.txt
rm ${outdir}/*_Sofk.txt
