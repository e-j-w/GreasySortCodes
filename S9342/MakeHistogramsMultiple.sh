#!/bin/bash

ADIR=AnalysisTrees
SMOLDIR=SmolFiles
HISTDIR=HistFiles
SORTCODE=SortDiagnosticsSMOL
FIRSTSUBRUN=0
LASTSUBRUN=10000000

if [ ! -d $HISTDIR ]; then
mkdir $HISTDIR
fi

if [ ! -d $SMOLDIR ]; then
mkdir $SMOLDIR
fi

if [ ! -d $SMOLDIR/reAB ]; then
mkdir $SMOLDIR/reAB
fi

if [ "$#" -eq 2 ]
then
	FIRSTRUN=$1
	LASTRUN=$2
	for i in $(seq $FIRSTRUN $LASTRUN);
	do
		echo "Making histograms for run $i" 
		./MakeHistograms.sh $i
	done

else
	echo "Provide the first and last run numbers as arguments!"
fi
