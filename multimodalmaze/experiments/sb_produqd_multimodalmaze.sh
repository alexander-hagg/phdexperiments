#!/bin/bash
#SBATCH --partition=any          # partition (queue) hpc, any, hpc3
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=4     # number of cores per node  
#SBATCH --mem=23G               # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=48:00:00          # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error=slurm.%j.err     # filename for STDERR
#SBATCH --export=ALL

# Run experiment
export JOBTMPDIR="/tmp/"$RANDOM
mkdir -p $JOBTMPDIR
rm -r $JOBTMPDIR/*
cd /home/ahagg2s/multimodalmaze/
matlab -nodisplay -nosplash -nodesktop -r "run_produqd"

