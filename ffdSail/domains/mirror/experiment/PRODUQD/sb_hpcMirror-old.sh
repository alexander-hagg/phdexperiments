#!/bin/sh
#PBS -N PRODUQD_Mirror
#PBS -q hpc1
#PBS -l nodes=1:ppn=16
#PBS -l walltime=72:00:00
#PBS -l vmem=100gb

# Necessary to pseudo-revert to old memory allocation behaviour
export MALLOC_ARENA_MAX=4

# Load Java, needed for parallel computing toolbox
# java/7, java/8 no noticable difference in terms of stability
module load java/default
module load cuda/default
module load matlab/R2017a

# Run experiment
cd $PBS_O_HOME/prodigi/ffdSail/domains/mirror/experiment/PRODIGI/
matlab -nodisplay -nosplash -nodesktop -r "mirror_hpcProduqd"
