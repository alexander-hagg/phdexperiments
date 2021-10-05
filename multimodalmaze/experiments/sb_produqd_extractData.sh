#!/bin/bash
#SBATCH --partition=any          # partition (queue) hpc, any, hpc3
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=32     # number of cores per node  
#SBATCH --mem=160G               # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=2:20:00          # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error=slurm.%j.err     # filename for STDERR
#SBATCH --export=ALL
#SBATCH --exclusive

# Run experiment
export JOBTMPDIR="/tmp/getData"
mkdir -p $JOBTMPDIR
cd /home/ahagg2s/multimodalmaze/
#matlab -nodisplay -nosplash -nodesktop -r "dataExtract('/scratch/ahagg2s/PROPHESAI/NH2', 2, 3, [30 90 150 210 270 330], 3, [1 2 3 9999],[25])"

matlab -nodisplay -nosplash -nodesktop -r "dataExtract('"$OUTPUTDIR"', "$NUMHIDDEN", "$NITERS", [30 90 150 210 270 330], 3, [1 2 3 9999],[25])"
# function dataExtract(folder,numHidden,numIterations,thetaValues,selCriterion,selValues,penalties)
