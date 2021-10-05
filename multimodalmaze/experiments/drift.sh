#!/bin/bash
echo "Running single ring experiment for GECCO 2019"
# Default parameters
export PBS_O_HOME=$home
export CFG_NCORES=30
export CFG_DATAPATH=/scratch/ahagg2s/GECCO2019/
export CFG_LOADFILE=0
export CFG_INPUTFILE=0
export CFG_RNN=1
export CFG_NUMHIDDEN=5
export CFG_SELCRITID=9999  # Default: no selection
export CFG_NITERS=1
export CFG_REPRESENTATION=planner
    
echo "Running experiment, selecting exits in ring 1 with representation"
echo $CFG_REPRESENTATION
export CFG_LOADFILE=1
export CFG_STARTITER=2
export CFG_NITERS=2
export CFG_INPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/baselines/"

# Experimental parameters
export CFG_SELECTIONMETHOD="individual" # class or individual
export CFG_DIMREDUCTION=1 # Use dimensionality reduction when comparing in similarity space
export CFG_ADJUSTOBJECTIVE=1 # Adjust objective function to include constraints
export CFG_SEEDING=0 # Set seeds. 

export ringID=1 # Ring 1, 2 or 3

# PLANNER export CFG_PWEIGHT=10
# CONTROLLER export CFG_PWEIGHT=200
if [ "CFG_REPRESENTATION" == "planner" ]; then
    export CFG_PWEIGHT=10
elif [ "CFG_REPRESENTATION" == "controller" ]; then
    export CFG_PWEIGHT=200
fi

mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION"
export CFG_SELCRITID=$(($ringID+1))

for experiment in 1 2 3; do
    if [ "$experiment" == "1" ]; then
        export TYPE=UDHM
        export CFG_ADJUSTOBJECTIVE=1
        export CFG_SEEDING=0 
    elif [ "$experiment" == "2" ]; then
        export TYPE=SEED
        export CFG_ADJUSTOBJECTIVE=0
        export CFG_SEEDING=1
    else
        export TYPE=COMB
        export CFG_ADJUSTOBJECTIVE=1
        export CFG_SEEDING=1
    fi

    # for mut in 0.000001 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.015 0.02 0.05 0.08 0.10 0.15 0.20; do
    # for mut in 0.00001 0.0001 0.001 0.10 0.12 0.15 0.20; do
    # PLANNER
    # for mut in 0.0001 0.0005 0.001 0.005 0.01 0.015 0.02 0.05 0.08; do
        export CFG_MUTRATE=$mut
        for exit1 in 1 2 3; do
             for a in 30 90 150 210 270 330; do
                export CFG_THETA="$a"
                export CFG_INPUTFILE="$CFG_THETA"
                echo $CFG_INPUTPATH$CFG_INPUTFILE 
                mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/""$TYPE"
                mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/""$TYPE""/exit""$exit1"
                mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/""$TYPE""/exit""$exit1""/""$CFG_MUTRATE"
                export CFG_OUTPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/""$TYPE""/exit""$exit1""/""$CFG_MUTRATE""/"
                export CFG_SELVAL=$exit1
                echo $CFG_SELVAL 
                if [ ! -f "$CFG_OUTPUTPATH""$a"".mat" ]; then
                    echo "$CFG_OUTPUTPATH""$a"".mat"" does not exist"
                    sbatch sb_produqd_multimodalmaze.sh
                else
                    echo "$CFG_OUTPUTPATH""$a"".mat"" does exist"
                fi
            done
        # sleep 3h
        done
        # sleep 25m
    done
done

# sbatch sb_produqd_extractDataMultiring.sh
