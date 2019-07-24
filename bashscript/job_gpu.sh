#!/bin/bash
#Submit this script with: sbatch thefilename

#asking for resources

#SBATCH --time=00:02:00   # walltime
#SBATCH --mem=20000M               # memory (per node)
#SBATCH --nodes=1

#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
##SBATCH --gpu-per-node=4
#SBATCH --sockets-per-node=1
#SBATCH --cores-per-socket=1
#SBATCH --ntasks-per-node=1 #for distributed memory mpi
#SBATCH --ntasks-per-socket=11
##SBATCH --cpus-per-gpu=1
#SBATCH --cpus-per-task=1  #for shared memory openmp

#SBATCH -J "speedtest"   # job name
#SBATCH --mail-user=hllin@caltech.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --qos=debug
#SBATCH --gres=gpu:1

SBATCH_GPUS_PER_NODE_local=1

# export environment
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
#export OMP_PROC_BIND=true
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# RETURN ENV
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
echo "SLURM_NTASKS_PER_NODE="$SLURM_NTASKS_PER_NODE
echo "SLURM_NTASKS_PER_SOCKET="$SLURM_NTASKS_PER_SOCKET
echo "SLURM_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "SLURM_GPUS_PER_NODE="$SLURM_GPUS_PER_NODE
echo "SBATCH_GRES="$SBATCH_GRES
echo "working directory = "$SLURM_SUBMIT_DIR
nvidia-smi
env
echo $PSM2_CUDA
lscpu
service nv_peer_mem status
lsmod | grep nv_peer_mem
lsmod | grep gdrdrv
# check CUDA_LAUNCH_BLOCKING dont set to 1
echo "CUDA_LAUNCH_BLOCKING="${CUDA_LAUNCH_BLOCKING} 

# ulimit -s 10240

#srun
SRUN_basic="srun --mpi=pmi2"
SRUNBIND_c="--cpu-bind=cores"
SRUNBIND_s="--cpu-bind=sockets"
#mpirun
OMPIRUN_basic="mpirun"
OMPIRUNBIND_c="--bind-to core"   #when the number of processes is <= 2
OMPIRUNBIND_s="--bind-to socket --map-by socket"  #when the number of processes is > 2

#====ompi environment====
#echo $OMPI_MCA_pml
#OMPI_MCA_pml=ucx
#export OMPI_MCA_pml

# ompi set gpu direct
export PSM2_GPUDIRECT=1
export OMPI_MCA_btl_openib_want_cuda_gdr=1
#ompi test
export OMPI_MCA_btl=^openib
#test
#export OMPI_MCA_btl="self,sm,tcp,^openib"

#====mvapich environment====
export MV2_USE_CUDA=1
export MV2_USE_GPUDIRECT_RDMA=1

#===control LAMMPS=============
LMP_CMD_kno="-k off"
LMP_CMD_kkgpu="-k on g ${SBATCH_GPUS_PER_NODE_local} -sf kk"
LMP_CMD_kkomp="-k on t ${SLURM_CPUS_PER_TASK} -sf kk"
LMP_CMD_kkgpuomp="-k on g ${SBATCH_GPUS_PER_NODE_local} t ${SLURM_CPUS_PER_TASK} -sf kk"
NEIGHDIRECT="-pk kokkos neigh half gpu/direct on"
NEIGHnoDIRECT="-pk kokkos neigh half gpu/direct off"

LMP_OMP="-pk omp ${SLURM_CPUS_PER_TASK} -sf omp"

#=======Script========================
LMP_INSCRIPT_sm="/home/hllin/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11"
LMP_INSCRIPT_kkbench_kokkos="/home/hllin/python/simulation/lammps_input_script/kokkos_bench/in.chute_kokkos"
LMP_INSCRIPT_kkbench_kokkos_print1="/home/hllin/python/simulation/lammps_input_script/kokkos_bench/in.chute_kokkos_print1atom"
LMP_INSCRIPT_ch="in.chute_print1atom"
LMP_INSCRIPT_lj="in.lj_print1atom"
LMP_INSCRIPT_myclinder="/home/hllin/python/simulation/lammps_input_script/in.lmpscript_20190114_v11"

LMP_INSCRIPT=${LMP_INSCRIPT_myclinder}

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

##===========ompi self using HPC ucx, must use gpu node============
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
${OMPIRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${OMPIRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kno} -in ${LMP_INSCRIPT}
${SRUN_basic} lmp ${LMP_CMD_kkgpu} ${NEIGHDIRECT} -in ${LMP_INSCRIPT}

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
##===========kkomp============
#LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load lmp/190723_master_ompi401_kkcpu ompi/4.0.1_yesucx_computenode

#check setting after load module
module list
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
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
#ompi setting check
ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
ompi_info --all | grep btl_openib_have_cuda_gdr
ompi_info --all | grep btl_openib_have_driver_gdr
#====run=====
${OMPIRUN_basic} ${OMPIRUNBIND_s} lmp ${LMP_CMD_OMP} -in ${LMP_INSCRIPT}
${SRUN_basic} ${SRUNBIND_s} lmp ${LMP_CMD_OMP} -in ${LMP_INSCRIPT}
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


#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half gpu/direct off -k on g 1 -sf kk -in ~/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in ~/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half gpu/direct off -k on g 1 -sf kk -in ~/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in ~/python/simulation/lammps_input_script/in.lmpscript_simple_20190114_v11 #${SLURM_CPUS_PER_TASK}

#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -pk kokkos neigh half -k on g 1 -sf kk -in in.chute #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda_debug/lmp -k on g 2 -pk kokkos neigh half -sf kk -in in.lj #${SLURM_CPUS_PER_TASK}
#srun --mpi=pmi2 ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda_debug/lmp -k off -in in.lj #${SLURM_CPUS_PER_TASK}
#mpirun ~/lmp_exe/20190712_mymodule_unstable_openmpi314_granular_kkcuda/lmp -k on g 1 -pk kokkos neigh half -sf kk -in in.lj #${SLURM_CPUS_PER_TASK}
'
echo "All Done!"
