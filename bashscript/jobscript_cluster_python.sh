#!/bin/bash
#Submit this script with: sbatch thefilename

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1
#SBATCH -J "20190529"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

##SBATCH --qos=debug

# RETURN ENV
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR

#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module list
#ulimit -s 10240
env

srun ~/simulation/script/script_at_dp.py

echo "All Done!"
