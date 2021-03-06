#!/bin/bash
#Submit this script with: sbatch thefilename

# test job
##SBATCH --qos=debug

# asking for resources
#SBATCH --time=20:30:00   # walltime
#SBATCH --mem-per-cpu=3G               # memory (per cpu)
#SBATCH --ntasks=32
#SBATCH --hint=nomultithread
#SBATCH --nodes=1
##SBATCH -B 2:16  # this also bound task to socket/core/thread since the slurm conf with affinity
##SBATCH --ntasks-per-node=28 #for distributed memory mpi

# ncpus per MPI task. choose ncpus processors per allocated CPU. (only use one)
##SBATCH --cpus-per-task=1  #for shared memory openmp

# job info
#SBATCH -J "H6"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

. ~/simulation/lammps_process/tools/set_environment/set_en_run_lmp_cpu.sh #--multi
