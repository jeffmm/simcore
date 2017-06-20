#!/bin/bash

#SBATCH --job-name=simcore
#SBATCH --time 00:30:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --output sc.out
#SBATCH --error sc.err
#SBATCH --account ucb-summit-smr
#SBATCH --qos=condo
#SBATCH --partition=shas

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=12
else
  omp_threads=1
fi
echo "starting sc job on $omp_threads threads"
export OMP_NUM_THREADS=$omp_threads

for i in `ls driving_mse2e_v*params.yaml`; do ./simcore -a $i; done

#time ./simcore params.yaml
