#!/bin/bash
#SBATCH --account CHENG-SL4-CPU #T2-CS061-CPU
#SBATCH --partition cclake #skylake
#SBATCH --nodes 1
#SBATCH --ntasks 8
#SBATCH --time 12:00:00

module purge
module load rhel7/default-peta4

source ~/.bashrc

mpirun -np 8 lmp_nnp_plumed < PREFIX-P-PRESSURE-T-TEMPERATURE.lmp > PREFIX-T-TEMPERATURE-P-PRESSURE.lmplog  &

wait

