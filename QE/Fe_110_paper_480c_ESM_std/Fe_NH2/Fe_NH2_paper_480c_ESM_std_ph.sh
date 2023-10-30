#!/bin/bash
#SBATCH --job-name=Fe_NH2_paper_480c_ESM_std_ph
#SBATCH --account=ucb316
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=249325M 
#SBATCH -t 48:00:000
#SBATCH --output=Fe_NH2_paper_480c_ESM_std_ph_%j.out
#SBATCH --error=Fe_NH2_paper_480c_ESM_std_ph_%j.err
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

srun --mpi=pmi2 -n 256 ph.x -npool 8 < Fe_NH2_paper_480c_ESM_std_ph.phi > Fe_NH2_paper_480c_ESM_std_ph.out