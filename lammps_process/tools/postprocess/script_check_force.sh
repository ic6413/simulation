#!/bin/bash

##### Functions
combine_all () {

    local START=$1
    local END=$2
    local STEPDIFF=$3

    # delete headers
    for i in $(eval echo "{$START..$END..$STEPDIFF}")
    do
    sed "-e 1,9d; s/^/$i /" output/single_all/dump.all.single.$i > output/single_all/tmp$i
    done

    cat $(eval echo "output/single_all/tmp{$START..$END..$STEPDIFF}") > output/single_all/dump.all.single.allstep
    
    for i in $(eval echo "{$START..$END..$STEPDIFF}")
    do
    rm -r output/single_all/tmp$i
    done


}

combine_trace () {

    local START=$1
    local END=$2
    local STEPDIFF=$3

    # delete headers
    for i in $(eval echo "{$START..$END..$STEPDIFF}")
    do
    sed "-e 1,9d; s/^/$i /" output/single_trace/dump.trace.single.$i > output/single_trace/tmp$i
    done

    cat $(eval echo "output/single_trace/tmp{$START..$END..$STEPDIFF}") > output/single_trace/dump.trace.single.allstep


    for i in $(eval echo "{$START..$END..$STEPDIFF}")
    do
    rm -r output/single_trace/tmp$i
    done


}
# check force calculation for atom i from step1 to step2

# check if the info of neighborhood of atom i from step1 to step2 have been output from LAMMPS
id_i=$1
step1=$2
step2=$3

use_trace_or_all=`python ~/simulation/lammps_process/python/script/script_check_dump_single_trace.py "$id_i" "$step1" "$step2"`
if [ "$use_trace_or_all" == "trace" ]; then
    # combine all data file from step1 to step2+1
    combine_trace ${step1} ${step2} "1"
elif [ "$use_trace_or_all" == "all" ]; then
    combine_all ${step1} ${step2} "1"
else
    echo "not trace not all"
fi

python ~/simulation/lammps_process/python/script/bash_to_python_check_force.py "$id_i" "$step1" "$step2" "$use_trace_or_all"
