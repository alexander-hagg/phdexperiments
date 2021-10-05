#!/bin/bash
echo "Running single ring experiment for GECCO 2019"
# Default parameters
export PBS_O_HOME=$home
export CFG_NCORES=32

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
export CFG_ADJUSTOBJECTIVE=0 # Adjust objective function to include constraints
export CFG_SELECTIONMETHOD="individual" # class or individual
export CFG_DIMREDUCTION=1 # Use dimensionality reduction when comparing in similarity space
export CFG_SEEDING=1 # Set seeds. 

export ringID=1 # Ring 1, 2 or 3
export CFG_PWEIGHT=1


mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION"
export CFG_SELCRITID=$(($ringID+1))

# sleep 360m

for exit1 in 1 2 3 9999; do
    mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/exit""$exit1"
    for a in 30 90 150 210 270 330; do
        export CFG_THETA="$a"
        export CFG_INPUTFILE="$CFG_THETA"
        echo $CFG_INPUTPATH$CFG_INPUTFILE 
        export CFG_OUTPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/exit""$exit1""/"     
        export CFG_SELVAL=$exit1
        echo $CFG_SELVAL
        export JOBTMPDIR="/tmp/""$a"
        if [ ! -f "$CFG_OUTPUTPATH""$a"".mat" ]; then
            echo "$CFG_OUTPUTPATH""$a"".mat"" does not exist"
            sbatch sb_produqd_multimodalmaze.sh
        else
            echo "$CFG_OUTPUTPATH""$a"".mat"" does exist, SKIPPING"
        fi
    done
    # sleep 140m
done

sbatch sb_produqd_extractDataMultiring.sh
