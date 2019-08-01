#!/bin/bash
#Submit this script with: sbatch thefilename

# asking for resources
#SBATCH --time=30:00:00   # walltime
#SBATCH --mem=20000M               # memory (per node)
#SBATCH --nodes=1
#SBATCH -B 2:16  # this also bound task to socket/core/thread since the slurm conf with affinity
#SBATCH --ntasks-per-node=32 #for distributed memory mpi

# ncpus per MPI task. choose ncpus processors per allocated CPU. (only use one)
#SBATCH --cpus-per-task=1  #for shared memory openmp

# test job
##SBATCH --qos=debug

# job info
#SBATCH -J "0729"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

. ~/python/simulation/bashscript/job_nodecpu.sh
