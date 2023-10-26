#!/bin/bash
#SBATCH --job-name=test_zdp
#SBATCH --account=ucb316
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH --array=0-499
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ketong_shao@berkeley.edu
## Command(s) to run:

declare -xr LOCAL_SCRATCH_DIR="/scratch/${USER}/job_${SLURM_JOB_ID}"
declare -xr ARRAY_TASK_DIR="${LOCAL_SCRATCH_DIR}/${SLURM_JOB_NAME}-${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

echo $LOCAL_SCRATCH_DIR
echo $ARRAY_TASK_DIR
echo $SLURM_SUBMIT_DIR

REALY_ID="`expr ${SLURM_ARRAY_TASK_ID}`"

export SLURM_EXPORT_ENV=ALL
module purge
module load cpu/0.15.4
module load gcc/10.2.0
module load anaconda3/2020.11
module load slurm 

mkdir -p "${ARRAY_TASK_DIR}"
cd "${ARRAY_TASK_DIR}"

python "${SLURM_SUBMIT_DIR}/Kinetic_funcs_entropy_array.py" "${REALY_ID}" "${SLURM_SUBMIT_DIR}" "${ARRAY_TASK_DIR}"

#cd "${LOCAL_SCRATCH_DIR}"

#cp -rp "${ARRAY_TASK_DIR}/${REALY_ID}_datastore" "${SLURM_SUBMIT_DIR}"