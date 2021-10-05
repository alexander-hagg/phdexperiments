#!/bin/sh
#PBS -N OpenFOAM_caseRunner
#PBS -q hpc1
#PBS -l nodes=1:ppn=12
#PBS -l walltime=48:00:00
#PBS -l vmem=120gb

# Necessary to pseudo-revert to old memory allocation behaviour
export MALLOC_ARENA_MAX=4

# Run experiment
cd $path
./Allclean && ./caseRunner.sh
