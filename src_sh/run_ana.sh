#!/bin/bash
# Optimized code from Subil/Omar -  OLCFHELP-22622 [cades-help] Unable to run Gromacs utilities (Check email: June-03-2025)
# Instructions to build GROMACS are in this email
# Date: July-02-2025

#SBATCH -A chem
#SBATCH -p burst
#SBATCH -t 06:00:00
#SBATCH -N 1                    # 1 nodes
#SBATCH --ntasks-per-node=1     # 1 tasks/node
#SBATCH -c 8                    # 8 cores/task
#SBATCH --mem=0
#SBATCH -J gengro_Fe3_m_0.1
#SBATCH -o ana_mol_0.1/out.%J
#SBATCH -e ana_mol_0.1/err.%J

# Load modules
module reset
module load intel mkl

# Export gcc path and num_threads
export LD_LIBRARY_PATH=/sw/cades-open/gcc/12.2.0/lib64:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Initializing jobs
cd $SLURM_SUBMIT_DIR
echo "begin job to analyze gro files  @start time: ${date}"
echo $PWD

outdir='ana_mol_0.1'
mkdir -p ${outdir}
wait

./ana.o anainp.txt
wait

cp anainp.txt ${outdir}
mv iondiff_* ${outdir}/
mv countiondiff_* ${outdir}/
mv autocorrcion_* ${outdir}/
mv autocorrpol_* ${outdir}/
mv rdf_* ${outdir}/
mv ciontype_* ${outdir}/
mv iontype_* ${outdir}/
mv clusttime_* ${outdir}/
mv scnt.txt ${outdir}/
mv all_neigh.txt ${outdir}/
mv clust_* ${outdir}/
mv catanneigh_* ${outdir}/
mv log_* ${outdir}/
mv polyionlist* ${outdir}/
mv polyclust_* ${outdir}/

echo "End of run.."
