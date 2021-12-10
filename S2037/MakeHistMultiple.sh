#!/bin/bash

firstrun=54736
lastrun=80000

if [ $# -eq 0 ]; then 
  echo "usage $0 <first run> <last run>"
  exit
fi

if [ $# -eq 1 ]; then 
  firstrun=$1
  lastrun=$1
fi

if [ $# -eq 2 ]; then
  firstrun=$1
  lastrun=$2
fi

# make histograms on a run by run basis
for run in `seq $firstrun $lastrun` ; do
  ./MakeHistograms.sh $run

done
