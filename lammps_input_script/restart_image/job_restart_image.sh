#!/bin/bash
START=$1
STEPDIFF=$2
END=$3
allrstimage=$4


if [ "${allrstimage}" == 1 ]
then
    echo "create image for all restart file"
    FILES=./output/rst/*
    for f in ${FILES}
    do
        lmp -in ${HOME}/simulation/lammps_input_script/restart_image/in.lmpscript_image_restart -var restartfile ${f} -log log.restartimage
    done
 else
    echo "create image for restart file at step " $1 "to " $3 "d_step is " $2
    for step in $(eval echo "{$START..$END..$STEPDIFF}")
    do
        lmp -in ${HOME}/simulation/lammps_input_script/restart_image/in.lmpscript_image_restart -var restartfile ./output/rst/restart.mpiio.${step} -log log.restartimage
    done   
fi

