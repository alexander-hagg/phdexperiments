#!/bin/bash
startFile="start"
stopFile="stop"
maxNumberCases=999

while true
    do
    FOLDER=$(($GPUID + 1))
	CASEFOLDER=$WORKFOLDER"/"$FOLDER
	mkdir $CASEFOLDER
	cd $CASEFOLDER	
	while true
		do
        sleep 10s

		if [ -f "$stopFile" ]
			then
            		cd $WORKFOLDER
		        echo "$stopFile found. Moving to next case"
			break
		fi

		if [ -f "$startFile" ]
			then
			echo "$startFile found: starting lettuce case."
			# Run experiment
			python building.py
		else
			echo -n "Waiting for $startFile..."
		fi

	done

done
