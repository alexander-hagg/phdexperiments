#!/bin/bash
# CRITERION 3: Use Exit X in Ring 1
echo "Running multiring experiment for GECCO 2019"

# module load matlab/default

# Default parameters
export PBS_O_HOME=$home
export CFG_NCORES=32

export CFG_DATAPATH=/scratch/ahagg2s/GECCO2019/
export CFG_LOADFILE=0
export CFG_INPUTFILE=0

export CFG_THETA=-60
export CFG_RNN=1
export CFG_NUMHIDDEN=5

export CFG_SELCRITID=9999  # Default: no selection
export CFG_NITERS=3

# CHANGE THIS IF NECESSARY
export CFG_ADJUSTOBJECTIVE=1

export BASELINE=false
# GET BASELINE
if $BASELINE;then
    echo "Getting baseline"
    export CFG_STARTITER=1
    echo "Getting base line"
	for a in 30 90 150 210 270 330; do
	   export CFG_NITERS=1
	   export CFG_THETA=$a
	   export CFG_SELVAL=$a+9999 
       export CFG_OUTPUTPATH="$CFG_DATAPATH""runs/baselines/"
	   export JOBTMPDIR="/tmp/"$a
	   sbatch sb_produqd_multimodalmaze.sh    
	done
    echo "Baseline complete. From now on load files from baseline or last round"
else 
    export CFG_LOADFILE=1
    export CFG_NITERS=1
    export CFG_INPUTPATH="$CFG_DATAPATH""runs/baselines/"

    #echo "Running experiment, selecting exits in ring 1"
    #export ringID=1
    #export CFG_STARTITER=$ringID
    # mkdir "$CFG_DATAPATH""runs"
    #export CFG_SELCRITID=$(($ringID+1))
    #for exit1 in 1 2 3; do
    #    for a in 30 90 150 210 270 330; do
    #        export CFG_THETA="$a"
    #        export CFG_INPUTFILE="$CFG_THETA"
    #        echo $CFG_INPUTPATH$CFG_INPUTFILE 
    #        mkdir "$CFG_DATAPATH""runs/exit""$exit1"
    #        export CFG_OUTPUTPATH="$CFG_DATAPATH""runs/exit""$exit1""/"     
    #        export CFG_SELVAL=$exit1
    #        echo $CFG_SELVAL
    #        export JOBTMPDIR="/tmp/""$a"
    #        sbatch sb_produqd_multimodalmaze.sh
    #    done
    #    sleep 25m
    #done

    echo "Running experiment, selecting exits in ring 2"
    export ringID=2
    export CFG_STARTITER=$ringID
    export CFG_SELCRITID=$(($ringID+1))

    for exit1 in 1 2 3; do
        export CFG_INPUTPATH="$CFG_DATAPATH""runs/exit""$exit1""/"
        for exit2 in 1 2 3; do
            for a in 30 90 150 210 270 330; do
                export CFG_THETA="$a"
                export CFG_INPUTFILE="$CFG_THETA"
                echo $CFG_INPUTPATH $CFG_INPUTFILE 
                mkdir "$CFG_DATAPATH""runs/exit""$exit1""/exit""$exit2"
                export CFG_OUTPUTPATH="$CFG_DATAPATH""runs/exit""$exit1""/exit""$exit2""/"     
                export CFG_SELVAL=$exit2
                echo $CFG_SELVAL
                export JOBTMPDIR="/tmp/""$a"
                sbatch sb_produqd_multimodalmaze.sh
            done
            sleep 25m
        done
    done
fi


# mkdir $OUTPUTDIR
# mv /scratch/ahagg2s/PROPHESAI/out* $OUTPUTDIR/
# sbatch sb_produqd_extractData.sh
