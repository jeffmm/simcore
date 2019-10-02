#!/bin/bash

#SBATCH --job-name=simcore
#SBATCH --time 00:30:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --output sc.out
#SBATCH --error sc.err
#SBATCH --account ucb-summit-smr
#SBATCH --qos=condo
#SBATCH --partition=shas
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jemo9179@colorado.edu

module purge
module load singularity/3.3.0

echo "Executing simcore job ${SLURM_JOB_ID} on ${SLURM_CPUS_PER_TASK} threads"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

singularity run simcore_latest.sif ./simcore.exe params.yaml
