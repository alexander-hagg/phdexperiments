#!/bin/bash
echo "Continuously try running new jobs"
for i in `seq 1 10`; do
    echo "Waiting to rerun"
    sleep 1h
    bash drift.sh
done

