#!/bin/bash
#SBATCH --job-name=N_480c
#SBATCH --account=fc_mllam
#SBATCH --partition=savio
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --output=N_480c_%j.out
#SBATCH --error=N_480c_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ketong_shao@berkeley.edu,jwang44@nd.edu
## Command(s) to run:

module load quantumespresso/6.4.1-qmcpack
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun pw.x -npool 4 < N_480c.pwi -in > N_480c.out