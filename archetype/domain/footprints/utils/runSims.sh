#!/bin/bash
conda activate lettuce
export GPUID=0

for i in {1..9}
do
   echo "Running simulation no. $i"
   cd $i
   python buildingDMD.py
   python DMD_multiple.py
   cd ..
done
