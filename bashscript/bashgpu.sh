#!/bin/bash
#Submit this script with: sbatch thefilename

# asking for resources
#SBATCH --time=00:03:00   # walltime
#SBATCH --mem=20000M               # memory (per node)
#SBATCH --nodes=1
##SBATCH 1:1:1  # this also bound task to socket/core/thread since the slurm conf with affinity
#SBATCH --ntasks-per-node=4 #for distributed memory mpi

# ncpus per MPI task, choose ncpus processors per allocated GPU or CPU. (only use one)
#SBATCH --cpus-per-task=1  #for shared memory openmp
# test job
#SBATCH --qos=debug

# job info
#SBATCH -J "speedtest"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# gpu per node
#SBATCH --gres=gpu:4
SBATCH_GPUS_PER_NODE_local=4

. ~/python/simulation/bashscript/job_nodegpu.sh
