#!/bin/bash
module load matlab/default

# Launch MAP-ELITES
echo 'MAP-Elites running on multimodal maze'
export PBS_O_HOME=$home
export CFG_NCORES=32

export CFG_DATAPATH=/scratch/ahagg2s/PROPHESAI/
export CFG_LOADFILE=0
export CFG_INPUTFILE=0
export CFG_SELECTIONMETHOD="individual" # class or individual
export CFG_DIMREDUCTION=0 # Use dimensionality reduction when comparing in similarity space
export CFG_ADJUSTOBJECTIVE=1 # Adjust objective function to include constraints
export CFG_SEEDING=0 # Set seeds. 

export CFG_RNN=0
export CFG_SELCRITID=9999  # Default: no selection
export CFG_NITERS=1

export CFG_REPRESENTATION=controller
export CFG_MUTRATE=0.1

export CFG_NUMHIDDEN=2
for a in 18 54 90 126 162 198 234 270 306; do
#  for a in 0 36 72 108 144 180 216 252 288 324; do
    export CFG_THETA=$a
    export SELVAL=$a+999999
    export CFG_OUTPUTPATH="$CFG_DATAPATH""PHD_2hidden/"
    sbatch sb_produqd_multimodalmaze.sh    
done

sleep 130m

export CFG_NUMHIDDEN=5
for a in 18 54 90 126 162 198 234 270 306; do
# for a in 0 36 72 108 144 180 216 252 288 324; do
    export CFG_THETA=$a
    export SELVAL=$a+999999
    export CFG_OUTPUTPATH="$CFG_DATAPATH""PHD_5hidden/"   
    sbatch sb_produqd_multimodalmaze.sh    
done
