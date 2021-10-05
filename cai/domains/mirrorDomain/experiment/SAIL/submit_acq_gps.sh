#!/bin/sh
#PBS -q default
#PBS -l nodes=1:ppn=12
#PBS -l walltime=10:00:00
#PBS -l vmem=60gb

export MALLOC_ARENA_MAX=4

module load java/7
module load matlab/R2017a

echo $PBS_JOBID
mkdir -p /scratch/$USER/tmp/$PBS_JOBID
mkdir -p /tmp/$PBS_JOBID

# Change to experiment folder
echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

ulimit -u 20000

# Run experiment
cd experiments/models
matlab -nodisplay -nosplash -nodesktop -r "acq_gps"

# Clean up temporary folder
rm -rf /tmp/$PBS_JOBID
