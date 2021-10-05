#!/bin/bash
#SBATCH --partition=gpu4        # partition (queue)
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=16     # number of cores per node
#SBATCH --mem=40G               # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=72:00:00          # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID) (OPTIONAL!!)
#SBATCH --error=slurm.%j.err     # filename for STDERR (OPTIONAL!!)
#SBATCH --export=ALL 		 # copy all environment variables from the user to the job (OPTIONAL!!)

module load matlab/R2019b

# Prepare Anaconda and Lettuce
source /home/$USER/anaconda3/bin/activate
conda init
conda activate lettuce

# Run experiment
cd /home/$USER/archetype/experiments

matlab -nodisplay -nosplash -nodesktop -r "runFinalSimulations"
