#!/bin/bash
#Submit this script with: sbatch thefilename

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1
##SBATCH --sockets-per-node=2
##SBATCH --cores-per-socket=16

#SBATCH -J "20190522_8"   # job name
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
module load openmpi/4.0.0  #mpich/3.2.1
module list
#ulimit -s 10240
env

. ./input_variable_script.sh

srun ~/work_script/lmp_mpi -var rst_from ${step1} -var runstep 1000000 -log log_${step1}.lammps -in in.lmpscript_20190114_v11

echo "All Done!"
