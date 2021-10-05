#!/bin/bash
for i in {1..1}
do
    echo /scratch/ahagg2s/$i
    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/$i submit_runSail.sh
done

