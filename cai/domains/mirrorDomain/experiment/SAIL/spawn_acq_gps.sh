#!/bin/bash
for i in {1..16}
do
    mkdir /scratch/ahagg2s/acq_parsec_2
    mkdir /scratch/ahagg2s/acq_parsec_2/$i
#    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/acquisition_vs_prediction_worst_vs_best_3/$i,RESTART=1 submit_runSail_acquisition_vs_prediction_worst_vs_best.sh
    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/acq_parsec_2/$i submit_acq_gps.sh
done

