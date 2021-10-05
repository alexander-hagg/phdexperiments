#!/bin/bash
#SBATCH --partition=any         # partition (queue) hpc, any, hpc3
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=8     # number of cores per node  
#SBATCH --mem=100G               # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=2:20:00           # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID)
#SBATCH --error=slurm.%j.err     # filename for STDERR
#SBATCH --export=ALL
#SBATCH --exclusive

# Run experiment
export JOBTMPDIR="/tmp/getData"
mkdir -p $JOBTMPDIR
cd /home/ahagg2s/multimodalmaze/

matlab -nodisplay -nosplash -nodesktop -r "analysisMutrate('/scratch/ahagg2s/GECCO2019/controller_MUTATION')"

