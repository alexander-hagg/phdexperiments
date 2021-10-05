#!/bin/sh
#PBS -q hpc
#PBS -l nodes=1:ppn=26
#PBS -l walltime=12:00:00
#PBS -l vmem=80gb

export MALLOC_ARENA_MAX=4

module load java/7
module load matlab/R2017a

echo $PBS_JOBID
mkdir -p /scratch/$USER/tmp/$PBS_JOBID
mkdir -p /tmp/$PBS_JOBID

# Change to experiment folder
cd $PBS_O_WORKDIR

ulimit -u 20000

# Run experiment
cd experiments/models
matlab -nodisplay -nosplash -nodesktop -r "analyzeRunSail"

# Clean up temporary folder
rm -rf /tmp/$PBS_JOBID
