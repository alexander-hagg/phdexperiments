#!/bin/bash
for i in {1..16}
do
    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/acq_parsec_2/$i submit_analyzeRunSail.sh
done


