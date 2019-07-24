#!/bin/bash
#Submit this script with: sbatch thefilename

# asking for resources
#SBATCH --time=00:02:00   # walltime
#SBATCH --mem=20000M               # memory (per node)
#SBATCH --nodes=1
#SBATCH -B=1:1:1
#SBATCH --ntasks-per-node=1 #for distributed memory mpi

# ncpus per MPI task. choose ncpus processors per allocated CPU. (only use one)
#SBATCH --cpus-per-task=1  #for shared memory openmp

# test job
#SBATCH --qos=debug

# job info
#SBATCH -J "speedtest"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# export environment
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_PROC_BIND=true
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
# ulimit -s 10240

# RETURN ENV
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
echo "SLURM_NTASKS_PER_NODE="$SLURM_NTASKS_PER_NODE
echo "SLURM_NTASKS_PER_SOCKET="$SLURM_NTASKS_PER_SOCKET
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "working directory = "$SLURM_SUBMIT_DIR
env
lscpu

#=================control MPI and LAMMPS job====================
#srun
SRUN_basic="srun --mpi=pmi2"
SRUNBIND_c="--cpu-bind=cores"
SRUNBIND_s="--cpu-bind=sockets"
#mpirun
OMPIRUN_basic="mpirun"
OMPIRUNBIND_c="--bind-to core"   #when the number of processes is <= 2
OMPIRUNBIND_s="--bind-to socket --map-by socket"  #when the number of processes is > 2

#====set ompi MCA environment====
export OMPI_MCA_btl=^openib
#test
#export OMPI_MCA_btl="self,sm,tcp,^openib"

#===control LAMMPS=============
LMP_CMD_kno="-k off"
LMP_CMD_kkomp="-k on t ${SLURM_CPUS_PER_TASK} -sf kk"
LMP_OMP="-pk omp ${SLURM_CPUS_PER_TASK} -sf omp"

#=======Script========================
LMP_INSCRIPT_sm="/home/hllin/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11"
LMP_INSCRIPT_kkbench_kokkos="/home/hllin/python/simulation/lammps_input_script/kokkos_bench/in.chute_kokkos"
LMP_INSCRIPT_kkbench_kokkos_print1="/home/hllin/python/simulation/lammps_input_script/kokkos_bench/in.chute_kokkos_print1atom"
LMP_INSCRIPT_ch="in.chute_print1atom"
LMP_INSCRIPT_lj="in.lj_print1atom"
LMP_INSCRIPT_myclinder="/home/hllin/python/simulation/lammps_input_script/in.lmpscript_20190114_v11"

LMP_INSCRIPT=${LMP_INSCRIPT_myclinder}

##===========kkomp============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load lmp/190723_master_ompi401_kkcpu ompi/4.0.1_yesucx_computenode

#check setting after load module
module list
#====run=====
${OMPIRUN_basic} ${OMPIRUNBIND_s} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} ${OMPIRUNBIND_s} lmp ${LMP_CMD_kkomp} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kkomp} -in ${LMP_INSCRIPT}


##===========user omp============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load lmp/190723_master_ompi401_useromp ompi/4.0.1_yesucx_computenode

#check setting after load module
module list
#====run=====
${OMPIRUN_basic} ${OMPIRUNBIND_s} lmp ${LMP_CMD_OMP} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_OMP} -in ${LMP_INSCRIPT}


echo "All Done!"
