#!/bin/bash
#SBATCH --partition=gpu4        # partition (queue)
#SBATCH --nodes=1                # number of nodes
#SBATCH --ntasks-per-node=12     # number of cores per node
#SBATCH --mem=30G               # memory per node in MB (different units with suffix K|M|G|T)
#SBATCH --time=72:00:00          # total runtime of job allocation (format D-HH:MM:SS; first parts optional)
#SBATCH --output=slurm.%j.out    # filename for STDOUT (%N: nodename, %j: job-ID) (OPTIONAL!!)
#SBATCH --error=slurm.%j.err     # filename for STDERR (OPTIONAL!!)
#SBATCH --export=ALL 		 # copy all environment variables from the user to the job (OPTIONAL!!)

module load cuda/default
echo $WORKFOLDER
echo $GPUID

# Prepare Anaconda and Lettuce
source /home/$USER/anaconda3/bin/activate
conda init
conda activate lettuce

bash /home/$USER/archetype/domain/footprints/exampleCallLettuce/lettuceWorker.sh
