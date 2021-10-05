#!/bin/sh
#PBS -N OpenFOAM_GroundTruth
#PBS -q hpc2
#PBS -l nodes=1:ppn=24
#PBS -l walltime=72:00:00
#PBS -l vmem=120gb

#---------------
# syntax: qsub -v fname=\'~/Code/data/ffdSail/sailFFD\',dname=\'~/Code/ffdSail/domains/foilFFD/\',dummy=1 getGroundTruth.sh
#
# or edit and run \'launch_getGroundTruth.sh\'
#
#---------------

# Necessary to pseudo-revert to old memory allocation behaviour
export MALLOC_ARENA_MAX=4

# Load Java, needed for parallel computing toolbox
# java/7, java/8 no noticable difference in terms of stability
module load java/default
module load cuda/default

# Run experiment
cd /home/agaier2m/Code/ffdSail/labBook/
matlab -nodisplay -nosplash -nodesktop -r "getGroundTruth($fname, $dname, $dummy)"
