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
export CFG_SELECTIONMETHOD="individual" # class or individual
export CFG_DIMREDUCTION=0 # Use dimensionality reduction when comparing in similarity space
export CFG_ADJUSTOBJECTIVE=1 # Adjust objective function to include constraints
export CFG_SEEDING=0 # Set seeds. 

# CHANGE THIS IF NECESSARY
export CFG_REPRESENTATION=controller
export CFG_MUTRATE=0.1

# Run experiment
echo "Getting baseline for representation "$CFG_REPRESENTATION
export CFG_STARTITER=1
mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION"
mkdir "$CFG_DATAPATH""$CFG_REPRESENTATION""/baselines"

for a in 30 90 150 210 270 330; do
   export CFG_THETA=$a
   export CFG_SELVAL=9999 
   export CFG_OUTPUTPATH="$CFG_DATAPATH""$CFG_REPRESENTATION""/baselines/"   
   sbatch sb_produqd_multimodalmaze.sh    
done
echo "Baseline complete."
