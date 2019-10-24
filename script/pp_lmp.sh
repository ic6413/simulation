#!/bin/bash
##### check number of input

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

##### Functions
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
            lmp -in ${HOME}/simulation/lmp_input/dump_restart/in.lmp_dump_image_at_restart -var restartfile ${f} -log log.restartimage
        done
    else
        echo "create image for restart file at step " $START "to " $END "d_step is " $STEPDIFF
        for step in $(eval echo "{$START..$END..$STEPDIFF}")
        do
            lmp -in ${HOME}/simulation/lmp_input/dump_restart/in.lmp_dump_image_at_restart -var restartfile ./output/rst/restart.mpiio.${step} -log log.restartimage
        done   
    fi
}

restart_movie () {
    
    local START_MOVIE=$1
    local STEPDIFF_MOVIE=$2
    local END_MOVIE=$3
    local n_framerate=$4

    rm -f /tmp/img*
    
    n_digit=$(bc -l <<< "l((($END_MOVIE - $START_MOVIE) / $STEPDIFF_MOVIE)+1)/l(10) + 1")
    n_digit=${n_digit%.*}

    x=1
    for i in $(eval echo "{$START_MOVIE..$END_MOVIE..$STEPDIFF_MOVIE}")
    do counter=$(printf %0"$n_digit"d $x)
        ln output/image/rst_all/all_image_"$i".jpg /tmp/img"$counter".jpg
        x=$(($x+1))
    done

    ffmpeg -f image2 -framerate ${n_framerate} -i /tmp/img%0${n_digit}d.jpg output/movie/movie_${START_MOVIE}_${END_MOVIE}.mov

    rm -f /tmp/img*

}

velocity_field () {
    
    local if_plot_to_last=$1
    local step1=$2
    local step2=$3
    local n_ave=$4

    if [ "$if_plot_to_last" == "1" ]; then
        echo "create velocity field at all step, n_ave is " $n_ave
    else
        echo "create velocity field at step " $step1 "to " $step2 "n_ave is " $n_ave
    fi
    ~/simulation/python/script/python_to_bash/script_velocity_field.py ${if_plot_to_last} ${step1} ${step2} ${n_ave}
}

wall_force () {
    
    local if_plot_to_last_wall=$1
    local step1_wall=$2
    local step2_wall=$3
    local n_ave_wall=$4

    if [ "$if_plot_to_last_wall" == "1" ]; then
        echo "create wall force at all step, n_ave_wall is " $n_ave_wall
    else
        echo "create wall force at step " $step1_wall "to " $step2 "n_ave_wall is " $n_ave_wall
    fi
    ~/simulation/python/script/python_to_bash/script_plot_wall_force.py ${if_plot_to_last_wall} ${step1_wall} ${step2_wall} ${n_ave_wall}
}

##### variable
interactive=
filename=~/sysinfo_page.html
START=
STEPDIFF=
END=
allrstimage=
START_MOVIE=
STEPDIFF_MOVIE=
END_MOVIE=
n_framerate=
if_plot_to_last=
step1=
step2=
n_ave=

###### command line option
while [ "$1" != "" ]; do
    case $1 in
        -f | --file )           shift
                                filename=$1
                                ;;
        -i | --interactive )    interactive=1
                                ;;
        --r2image )             if_r2image=1
                                shift
                                START=$1
                                shift
                                STEPDIFF=$1
                                shift
                                END=$1
                                shift
                                allrstimage=$1
                                ;;
        --r2movie )             if_r2movie=1
                                shift
                                START_MOVIE=$1
                                shift
                                STEPDIFF_MOVIE=$1
                                shift
                                END_MOVIE=$1
                                shift
                                n_framerate=$1
                                ;;
        --vfield )              if_vfield=1
                                shift
                                if_plot_to_last=$1
                                shift
                                step1=$1
                                shift
                                step2=$1
                                shift
                                n_ave=$1
                                ;;
        --wall )                if_wall=1
                                shift
                                if_plot_to_last_wall=$1
                                shift
                                step1_wall=$1
                                shift
                                step2_wall=$1
                                shift
                                n_ave_wall=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

###### Main

# check file name
echo "output file = $filename"

# if interactive then input variable
if [ "$interactive" == "1" ]; then
    echo "interactive is on"
    
    echo "Type the if_r2image that you want (1 or 0), followed by [ENTER]:"
    read if_r2image
    if [ "$if_r2image" == "1" ]; then
        echo "Type the allrstimage that you want (1 or 0), followed by [ENTER]:"
        read allrstimage
        echo "Type the START that you want, followed by [ENTER]:"
        read START
        echo "Type the STEPDIFF that you want, followed by [ENTER]:"
        read STEPDIFF
        echo "Type the END that you want, followed by [ENTER]:"
        read END
    fi

    echo "Type the if_r2movie that you want (1 or 0), followed by [ENTER]:"
    read if_r2movie
    if [ "$if_r2movie" == "1" ]; then
        echo "Type the n_framerate that you want, followed by [ENTER]:"
        read n_framerate
        echo "Type the START_MOVIE that you want, followed by [ENTER]:"
        read START_MOVIE
        echo "Type the STEPDIFF_MOVIE that you want, followed by [ENTER]:"
        read STEPDIFF_MOVIE
        echo "Type the END_MOVIE that you want, followed by [ENTER]:"
        read END_MOVIE
    fi

    echo "Type the if_vfield that you want (1 or 0), followed by [ENTER]:"
    read if_vfield
    if [ "$if_vfield" == "1" ]; then
        echo "Type the if_plot_to_last that you want (1 or 0), followed by [ENTER]:"
        read if_plot_to_last
        echo "Type the step1 that you want, followed by [ENTER]:"
        read step1
        echo "Type the step2 that you want, followed by [ENTER]:"
        read step2
        echo "Type the n_ave that you want, followed by [ENTER]:"
        read n_ave
    fi

    echo "Type the if_wall that you want (1 or 0), followed by [ENTER]:"
    read if_wall
    if [ "$if_wall" == "1" ]; then
        echo "Type the if_plot_to_last_wall that you want (1 or 0), followed by [ENTER]:"
        read if_plot_to_last_wall
        echo "Type the step1_wall that you want, followed by [ENTER]:"
        read step1_wall
        echo "Type the step2_wall that you want, followed by [ENTER]:"
        read step2_wall
        echo "Type the n_ave_wall that you want, followed by [ENTER]:"
        read n_ave_wall
    fi

else
    echo "interactive off"
fi

# r2image start
if [ "$if_r2image" == "1" ]; then
	echo "if_r2image is on"
    restart_image ${START} ${STEPDIFF} ${END} ${allrstimage}
else
	echo "if_r2image is off"
fi

# r2movie start
if [ "$if_r2movie" == "1" ]; then
	echo "if_r2movie is on"
    restart_movie ${START_MOVIE} ${STEPDIFF_MOVIE} ${END_MOVIE} ${n_framerate}
else
	echo "if_r2movie is off"
fi

# vfield start
if [ "$if_vfield" == "1" ]; then
	echo "if_vfield is on"
    velocity_field ${if_plot_to_last} ${step1} ${step2} ${n_ave}
else
	echo "if_vfield is off"
fi

# wall start
if [ "$if_wall" == "1" ]; then
	echo "if_wall is on"
    wall_force ${if_plot_to_last_wall} ${step1_wall} ${step2_wall} ${n_ave_wall}
else
	echo "if_wall is off"
fi