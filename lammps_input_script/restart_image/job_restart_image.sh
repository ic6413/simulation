#!/bin/bash
START=0
STEPDIFF=1
END=2
allrstimage=0

if [ "${allrstimage}" == 1 ]
then
    FILES=./output/rst/*
    for f in ${FILES}
    do
        lmp -in ${HOME}/simulation/lammps_input_script/restart_image/in.lmpscript_image_restart -var restartfile ${f}
    done
 else
    for step in $(eval echo "{$START..$END..$STEPDIFF}")
    do
        lmp -in ${HOME}/simulation/lammps_input_script/restart_image/in.lmpscript_image_restart -var restartfile ./output/rst/restart.mpiio.${step}
    done   
fi

