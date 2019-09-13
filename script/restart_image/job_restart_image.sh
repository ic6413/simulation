#!/bin/bash


if [ "$1" != "" ]; then
    echo "Positional parameter 1 contains something"
else
    echo "Positional parameter 1 is empty"
fi

if [ $# -gt 0 ]; then
    echo "Your command line contains $# arguments"
else
    echo "Your command line contains no arguments"
fi

interactive=
filename=~/sysinfo_page.html
START=0
STEPDIFF=0
END=0
allrstimage=1

while [ "$1" != "" ]; do
    case $1 in
        -f | --file )           shift
                                filename=$1
                                ;;
        -i | --interactive )    interactive=1
                                ;;
        -r | --reimage )        shift
                                START=$1
                                shift
                                STEPDIFF=$1
                                shift
                                END=$1
                                shift
                                allrstimage=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# Test code to verify command line processing

if [ "$interactive" == "1" ]; then
	echo "interactive is on"
else
	echo "interactive is off"
fi
echo "output file = $filename"


restart_image () {
    
    local START=$1
    local STEPDIFF=$2
    local END=$3
    local allrstimage=$4


    if [ "$allrstimage" == "1" ]; then
        echo "create image for all restart file"
        FILES=./output/rst/*
        for f in ${FILES}
        do
            lmp -in ${HOME}/simulation/lammps_input_script/dump_restart/in.lmp_dump_image_at_restart -var restartfile ${f} -log log.restartimage
        done
    else
        echo "create image for restart file at step " $1 "to " $3 "d_step is " $2
        for step in $(eval echo "{$START..$END..$STEPDIFF}")
        do
            lmp -in ${HOME}/simulation/lammps_input_script/dump_restart/in.lmp_dump_image_at_restart -var restartfile ./output/rst/restart.mpiio.${step} -log log.restartimage
        done   
    fi
}
if [ "$interactive" == "1" ]; then
    echo "Type the allrstimage that you want (1 or 0), followed by [ENTER]:"
    read allrstimage
    echo "Type the START that you want, followed by [ENTER]:"
    read START
    echo "Type the STEPDIFF that you want, followed by [ENTER]:"
    read STEPDIFF
    echo "Type the END that you want, followed by [ENTER]:"
    read END
else
    echo "interactive off"
fi

restart_image ${START} ${STEPDIFF} ${END} ${allrstimage}
