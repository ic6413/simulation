#!/bin/bash
#Submit this script with: sbatch thefilename

#SBATCH --time=00:01:00   # walltime
##SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=1

#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH -J "20190708"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --qos=debug
# RETURN ENV
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR


#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load openmpi/2.1.2 #openmpi/3.0.0  mpich/3.2.1
module list
#ulimit -s 10240

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
env
srun --cpu_bind=verbose,sockets ~/lmp_exe/20190708_unstable_openmpi212_nokkcuda_mpi_granular/lmp -in in.chute_regular

echo "All Done!"
