#!/bin/bash

#qsub -v fname=\'~/Code/data/ffdSail/sailFFD10\',dname=\'~/Code/ffdSail/domains/foilFFD/\',dummy=0 getGroundTruth.sh
#qsub -v fname=\'~/Code/data/ffdSail/sailFFD20\',dname=\'~/Code/ffdSail/domains/foilFFD/\',dummy=0 getGroundTruth.sh
#qsub -v fname=\'~/Code/data/ffdSail/sailParsec10\',dname=\'~/Code/ffdSail/domains/airfoil/\',dummy=0 getGroundTruth.sh
#qsub -v fname=\'~/Code/data/ffdSail/sailParsec20\',dname=\'~/Code/ffdSail/domains/airfoil/\',dummy=0 getGroundTruth.sh
qsub -v fname=\'~/Code/data/ffdSail/sailCPPN01\',dname=\'~/Code/ffdSail/domains/foilFFD/\',dummy=1 getGroundTruth.sh