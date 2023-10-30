#!/bin/bash
#SBATCH --job-name=Fe_NH_paper_480c_ESM_std
#SBATCH --account=ucb316
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --mem=249325M 
#SBATCH -t 48:00:00
#SBATCH --output=Fe_NH_paper_480c_ESM_std_%j.out
#SBATCH --error=Fe_NH_paper_480c_ESM_std_%j.err
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

srun --mpi=pmi2 -n 256 pw.x -npool 8 < Fe_NH_paper_480c_ESM_std.pwi > Fe_NH_paper_480c_ESM_std.out