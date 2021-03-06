# export environment
#export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_PROC_BIND=true
#export OMP_PROC_BIND=spread
#export OMP_PLACES=threads
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
#export OMPI_MCA_btl=^openib
#test
#export OMPI_MCA_btl="self,sm,tcp,^openib"

#===control LAMMPS=============
LMP_CMD_kno="-k off"
LMP_CMD_kkomp="-k on t ${SLURM_CPUS_PER_TASK} -sf kk"
LMP_OMP="-pk omp ${SLURM_CPUS_PER_TASK} -sf omp"


LMP_INSCRIPT=./in.override


##===========normal do not bind socket============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load openmpi/4.0.1
module load lmp/20200402_stable3Mar2020_openmpi401_gran_mpiio_useromp #190809unstable_openmpi401_mpiio #20200402_stable3Mar2020_openmpi401_gran_mpiio_useromp # 190809unstable_openmpi401_mpiio  190916master_openmpi401_omp_mpiio  #190809unstable_openmpi401_mpiio
#check setting after load module
module list
${OMPIRUN_basic} lmp -in ${LMP_INSCRIPT}

'
##===========normal do not bind socket============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load ompi/4.0.1_yesucx_computenode
module load lmp/190729_master_selfompi401ucx_normal_mpiio 
#check setting after load module
module list
${OMPIRUN_basic} lmp -in ${LMP_INSCRIPT}
${OMPIRUN_basic} ${OMPIRUNBIND_c} lmp -in ${LMP_INSCRIPT}
'

echo "All Done!"


