#!/bin/bash
#SBATCH --job-name=H2_480c_ESM_VA8_ph
#SBATCH --account=ucb316
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=120
#SBATCH --cpus-per-task=1
#SBATCH --mem=240000M 
#SBATCH -t 48:00:00
#SBATCH --output=H2_480c_ESM_VA8_ph_%j.out
#SBATCH --error=H2_480c_ESM_VA8_ph_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ketong_shao@berkeley.edu,jwang44@nd.edu
## Command(s) to run:

export SLURM_EXPORT_ENV=ALL
module purge
module load cpu/0.15.4
module load gcc/9.2.0
module load openmpi/3.1.6
module load quantum-espresso/6.7.0-openblas
module load slurm

srun -n 120 ph.x -npool 8 < H2_480c_ESM_VA8_ph.phi > H2_480c_ESM_VA8_ph.out