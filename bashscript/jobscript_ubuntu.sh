#!/bin/bash
. ./input_variable_script.sh

mpirun -np 4 ~/lammps_simulation/work_script/lmp_mpi_ubuntu -var rst_from ${step1} -var runstep ${numbersteptorun} -log log_${step1}.lammps -in in.lmpscript_20190114_v11
