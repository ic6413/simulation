
# export environment
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
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

# RETURN env related to GPU
echo "SBATCH_GRES="$SBATCH_GRES
echo "PSM2_CUDA="$PSM2_CUDA
# check CUDA_LAUNCH_BLOCKING not set to 1
echo "CUDA_LAUNCH_BLOCKING="${CUDA_LAUNCH_BLOCKING}  
service nv_peer_mem status
lsmod | grep nv_peer_mem
lsmod | grep gdrdrv
nvidia-smi

#=================control MPI and LAMMPS job====================
#srun
SRUN_basic="srun --mpi=pmi2 --accel-bind=g"
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

#====set ompi GPU environment====
export PSM2_GPUDIRECT=1
export OMPI_MCA_btl_openib_want_cuda_gdr=1

#====set mvapich2 GPU environment====
export MV2_USE_CUDA=1
export MV2_USE_GPUDIRECT_RDMA=1

#===control LAMMPS=============
LMP_CMD_kno="-k off"
LMP_CMD_kkgpu="-k on g ${SBATCH_GPUS_PER_NODE_local} -sf kk"
LMP_CMD_kkomp="-k on t ${SLURM_CPUS_PER_TASK} -sf kk"
LMP_CMD_kkgpuomp="-k on g ${SBATCH_GPUS_PER_NODE_local} t ${SLURM_CPUS_PER_TASK} -sf kk"
NEIGHDIRECT="-pk kokkos neigh half gpu/direct on" #newton off comm host
NEIGHnoDIRECT="-pk kokkos neigh half gpu/direct off"
LMP_OMP="-pk omp ${SLURM_CPUS_PER_TASK} -sf omp"

#=======Script========================
LMP_INSCRIPT_sm="/home/hllin/simulation/lmp_input/in.lmpscript_simple_20190114_v11"
LMP_INSCRIPT_kkbench_kokkos="/home/hllin/simulation/lmp_input/kokkos_bench/in.chute_kokkos"
LMP_INSCRIPT_kkbench_kokkos_print1="/home/hllin/simulation/lmp_input/kokkos_bench/in.chute_kokkos_print1atom"
LMP_INSCRIPT_ch="in.chute_print1atom"
LMP_INSCRIPT_lj="in.lj_print1atom"
LMP_INSCRIPT_mycylinder="/home/hllin/simulation/lmp_input/in.lmpscript_cylinder"
LMP_INSCRIPT_infolder="in.override"
LMP_INSCRIPT=${LMP_INSCRIPT_infolder}

:'
##===========ompi self usually slow============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
#module load cuda/10.0 ompi/4.0.1_cuda10.0_yesucx lmp/190719master_ompi401_cuda10_yesucx
module load cuda/10.0 ompi/4.0.1_cuda10.0_noucx lmp/190718master_ompi401_cuda10_noucx
#check setting after load module
module list
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
#====run=====
${OMPIRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT} # 14coreMPI best
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHnoDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHnoDIRECT} -in ${LMP_INSCRIPT}
'

##===========kkcudaomp, ompi self using HPC ucx, must use gpu node============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load cuda/10.0 ucx/1.5.1_cuda-10.0 ompi/4.0.1_cuda10.0_yesucx lmp/190724_master_ucxHPC_selfompi401_kk_cuda10_omp

#check setting after load module
module list
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
#====run=====
:'
${OMPIRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHnoDIRECT} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
'
${SRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHnoDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHDIRECT} comm host -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHnoDIRECT} comm host -in ${LMP_INSCRIPT}

##===========kkcuda, ompi self using HPC ucx, must use gpu node============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load cuda/10.0 ucx/1.5.1_cuda-10.0 ompi/4.0.1_cuda10.0_yesucx lmp/190723_master_ucxHPC_selfompi401_cuda10

#check setting after load module
module list
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
#====run=====
:'
${OMPIRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpuomp} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
'
${SRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHnoDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} comm host -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHnoDIRECT} comm host -in ${LMP_INSCRIPT}

:'
##===========ompi HPC can not use srun============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load cuda/10.0 ucx/1.5.1_cuda-10.0 openmpi/4.0.1_cuda-10.0 lmp/190723_master_HPC_ompi401_cuda10_ucx151

#check setting after load module
module list
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
#====run=====
${OMPIRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
'


:'
##===========mvapich2============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load cuda/10.0 mvapich2/2.3.1-gdr-cuda_10.0 lmp/190719master_ompi401_cuda10_yesucx
#check setting after load module
module list
#====run=====
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kkomp} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_kkgpuomp} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
'

:'
module purge
module load cuda/8.0 openmpi/3.1.4
#module load cuda/10.0 ompi/3.1.4_cuda10.0 #openmpi/3.0.0  mpich/3.2.1
#module load ucx/1.6_cuda10.0_nonuma cuda/10.0 ompi/3.1.4_ucx1.6_cuda10.0_nonuma
module list
${MPIRUN_SCRIPT} ${LMP_EXE_clusteropenmpi} ${LMP_CMD_kno} ${LMP_INSCRIPT} 
${MPIRUN_SCRIPT} ${LMP_EXE_clusteropenmpi} ${LMP_CMD_kno} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_clusteropenmpi} ${LMP_CMD_kyes} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_clusteropenmpi} ${LMP_CMD_kyes_directon} ${LMP_INSCRIPT}

module purge
module load cuda/10.0 ompi/3.1.4_cuda10.0 #openmpi/3.0.0  mpich/3.2.1
module list
${MPIRUN_SCRIPT} ${LMP_EXE_cuda10} ${LMP_CMD_kno} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_cuda10} ${LMP_CMD_kno} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_cuda10} ${LMP_CMD_kyes} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_cuda10} ${LMP_CMD_kyes_directon} ${LMP_INSCRIPT}


module purge
module load ucx/1.6_cuda10.0_nonuma cuda/10.0 ompi/3.1.4_ucx1.6_cuda10.0_nonuma
module list
${MPIRUN_SCRIPT} ${LMP_EXE_ucx16_cuda100_nonuma} ${LMP_CMD_kno} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_ucx16_cuda100_nonuma} ${LMP_CMD_kno} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_ucx16_cuda100_nonuma} ${LMP_CMD_kyes} ${LMP_INSCRIPT}
${MPIRUN_SCRIPT} ${LMP_EXE_ucx16_cuda100_nonuma} ${LMP_CMD_kyes_directon} ${LMP_INSCRIPT}

#${SRUN_SCRIPT} ${LMP_EXE_cuda10} ${LMP_CMD_kyes_directon} ${LMP_INSCRIPT_kkbench_kokkos}
#${SRUN_SCRIPT} ${LMP_EXE_ucx16_cuda100_nonuma} ${LMP_CMD_kyes} ${LMP_INSCRIPT_kkbench_kokkos}


#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half gpu/direct off -k on g 1 -sf kk -in ~/simulation/lmp_input/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in ~/simulation/lmp_input/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half gpu/direct off -k on g 1 -sf kk -in ~/simulation/lmp_input/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in ~/simulation/lmp_input/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}

#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in in.chute #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda_debug/lmp -k on g 2 -pk kokkos neigh half -sf kk -in in.lj #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda_debug/lmp -k off -in in.lj #${SLURM_CPUS_PER_TASK}
#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -k on g 1 -pk kokkos neigh half -sf kk -in in.lj #${SLURM_CPUS_PER_TASK}
'
echo "All Done!"
