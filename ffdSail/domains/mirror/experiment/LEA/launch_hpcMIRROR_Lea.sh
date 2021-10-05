#!/bin/bash

# Create base OpenFOAM cases and launch case runners
 nHpc1Cases=3
 nHpc2Cases=0
 destFolderName="/scratch/ahagg2s/sailCFD/"

# HPC1 with 12 Cores (1 job per node)
baseFolderName="/home/ahagg2s/prodigi/ffdSail/domains/mirror/pe/ofTemplates/hpc1_mirror_only"
for (( i=1; i<=$nHpc1Cases; i++ ))
do
	caseName=$destFolderName"case$i"
	echo $caseName
	cp -TR $baseFolderName $caseName
	qsub -j oe -v path="$caseName" $caseName/submit.sh
done 

# HPC2 with 12 Cores (2 jobs per node)
baseFolderName="/home/ahagg2s/prodigi/ffdSail/domains/mirror/pe/ofTemplates/hpc2_mirror_only"
for (( i=nHpc1Cases+1; i<=$nHpc1Cases+$nHpc2Cases; i++ ))
do
	caseName=$destFolderName"case$i"
	echo $caseName
	cp -TR $baseFolderName $caseName
	qsub -j oe -v path="$caseName" $caseName/submit.sh
done 

# Launch SAIL
cases=$(($nHpc1Cases + $nHpc2Cases))
echo 'PRODUQD Main Script'
qsub -N PRODUQD_FFD_MIRROR -v encoding=\'ffd\',nCases=3,startCase=1 sb_hpcMIRROR_Lea.sh
