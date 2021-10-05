#!/bin/bash
echo "Parameterization of penalty weight for GECCO 2019"
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

export CFG_REPRESENTATION=controller
    
echo "Running experiment, selecting exits in ring 1 with representation"
echo $CFG_REPRESENTATION

export CFG_LOADFILE=1
export CFG_STARTITER=2
export CFG_NITERS=2
export CFG_INPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/baselines/"

# Experimental parameters
export ringID=1 # Ring 1, 2 or 3
export CFG_SELECTIONMETHOD="individual" # class or individual
export CFG_ADJUSTOBJECTIVE=1 # Adjust objective function to include constraints
export CFG_DIMREDUCTION=1 # Use dimensionality reduction when comparing in similarity space
export CFG_SEEDING=0 # Set seeds. 

export CFG_MUTRATE=0.01

mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION"
export CFG_SELCRITID=$(($ringID+1))
export CFG_SELVAL=1

for pweight in 0 0.5 1 2.5 5 10 25 50 75 100 150 200 250 500; do
    export CFG_PWEIGHT=$pweight
    for a in 30 90 150 210 270 330; do
        export CFG_THETA="$a"
        export CFG_INPUTFILE="$CFG_THETA"
        echo $CFG_INPUTPATH$CFG_INPUTFILE 
        mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/pweight""$pweight"
        export CFG_OUTPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/pweight""$pweight""/"
        if [ ! -f "$CFG_OUTPUTPATH""$a"".mat" ]; then
            echo "$CFG_OUTPUTPATH""$a"".mat"" does not exist"
            sbatch sb_produqd_multimodalmaze.sh
        else
            echo "$CFG_OUTPUTPATH""$a"".mat"" does exist"
        fi
    done
done

# sbatch sb_produqd_extractDataParameterization.sh
