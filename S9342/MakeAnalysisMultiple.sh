#!/bin/bash

ADIR=AnalysisTrees
SMOLDIR=SmolFiles
FIRSTSUBRUN=0
LASTSUBRUN=10000000

if [ "$#" -eq 2 ]
then
	FIRSTRUN=$1
	LASTRUN=$2
	for i in $(seq $FIRSTRUN $LASTRUN);
	do
		echo "Making analysis trees for run $i" 
		./MakeAnalysis.sh $i
	done

else
	echo "Provide the first and last run numbers as arguments!"
fi
