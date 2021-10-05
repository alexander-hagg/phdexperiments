#!/bin/bash
conda activate lettuce
export numGPUs=4
export WORKFOLDER=/scratch/$USER/ARCHETYPE/workfolderCLUSTER
mkdir $WORKFOLDER
cd "$WORKFOLDER" 

# Start caserunners in each directory
for ((i=0; i<numGPUs; i++)); do
	export GPUID=$i
	echo "Starting GPU/Worker: $GPUID"
	sbatch /home/$USER/archetype/domain/footprints/exampleCallLettuce/submitLettuceWorker.sh &
done

sbatch /home/$USER/archetype/domain/footprints/exampleCallLettuce/submitLettuceSAILMainScript.sh &
