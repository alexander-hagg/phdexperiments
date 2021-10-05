#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=24
#PBS -l walltime=48:00:00
#PBS -l vmem=120gb

export MALLOC_ARENA_MAX=4

echo $PBS_JOBID
mkdir -p /scratch/$USER/tmp/$PBS_JOBID
mkdir -p /tmp/$PBS_JOBID

# Change to experiment folder
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# Run experiment
./Allclean
./Allrun

# Clean up temporary folder
rm -rf /tmp/$PBS_JOBID
