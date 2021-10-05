#!/bin/bash

# Create base OpenFOAM cases and launch case runners
 nHpc1Cases=0
 nHpc2Cases=2
 destFolderName="/scratch/ahagg2s/sailCFD/"

# HPC1 with 12 Cores (1 job per node)
baseFolderName="/home/ahagg2s/sail/domains/mirror/pe/ofTemplates/hpc1/"
for (( i=1; i<=$nHpc1Cases; i++ ))
do
	caseName=$destFolderName"case$i"
	echo $caseName
	cp -TR $baseFolderName $caseName
	qsub -j oe -v path="$caseName" $caseName/submit.sh
done 

# HPC2 with 12 Cores (2 jobs per node)
baseFolderName="/home/ahagg2s/sail/domains/mirror/pe/ofTemplates/hpc2/"
for (( i=nHpc1Cases+1; i<=$nHpc1Cases+$nHpc2Cases; i++ ))
do
	caseName=$destFolderName"case$i"
	echo $caseName
	cp -TR $baseFolderName $caseName
	qsub -j oe -v path="$caseName" $caseName/submit.sh
done 

# Launch SAIL
cases=$(($nHpc1Cases + $nHpc2Cases))
echo 'Mirror Test Main Script'
qsub -N FOAM_TEST-mirror -v encoding=\'ffd\',nCases=2,startCase=1 /home/ahagg2s/sail/domains/mirror/experiment/sb_hpcMirrorFoamTest.sh
