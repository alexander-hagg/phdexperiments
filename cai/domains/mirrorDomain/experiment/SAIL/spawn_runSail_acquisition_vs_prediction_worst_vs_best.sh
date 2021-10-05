#!/bin/bash
for i in {1..1}
do
    mkdir /scratch/ahagg2s/acq_GP
#    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/acquisition_vs_prediction_worst_vs_best_3/$i,RESTART=1 submit_runSail_acquisition_vs_prediction_worst_vs_best.sh
    qsub -v EXPERIMENTPATH=/scratch/ahagg2s/acq_GP/$i submit_runSail_acquisition_vs_prediction_worst_vs_best.sh
done

