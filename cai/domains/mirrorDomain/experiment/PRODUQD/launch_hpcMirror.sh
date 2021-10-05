#!/bin/bash
of240
module load matlab/default

# Create base OpenFOAM cases and launch case runners
 nHpc1Cases=10
 nHpc2Cases=0
 destFolderName="/scratch/ahagg2s/sailCFD/"

# HPC1 with 32 Cores (1 job per node)
baseFolderName="/home/ahagg2s/produqd/ffdSail/domains/mirror/pe/ofTemplates/hpc1_mirror_only"
for (( i=1; i<=$nHpc1Cases; i++ ))
do
	caseName=$destFolderName"case$i"
	echo $caseName
	cp -TR $baseFolderName $caseName
	# qsub -j oe -v path="$caseName" $caseName/submit.sh
 	export PBS_WORKDIR="$caseName"
	sbatch $caseName/submit.sh
done 

# Launch SAIL
cases=$(($nHpc1Cases))
echo 'PRODUQD Main Script'
# qsub -N PRODUQD_FFD_MIRROR -v encoding=\'ffd\',nCases=nHpc1Cases,startCase=1 sb_hpcMirror.sh
export encoding=\'ffd\'
export nCases=nHpc1Cases
export startCase=1
export PBS_O_HOME=$home
sbatch sb_hpcMirror.sh
