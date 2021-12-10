#!/bin/bash

declare -a FILE_LIST # array to hold file names
HIST_DIR=./HistFiles

firstrun=55775
lastrun=80000

if [ $# -eq 0 ]; then 
  echo "Sums all run histograms into single file"
  echo "usage $0 <first run> <last run>"
  exit
fi

if [ $# -eq 1 ]; then 
  echo "Sums all run histograms into single file"
  echo "usage $0 <first run> <last run>"
  exit
fi

if [ $# -eq 2 ]; then
  firstrun=$1
  lastrun=$2
fi

# summing all runs
for run in `seq $firstrun $lastrun` ; do
  if [[ -e $HIST_DIR/Runs_"$firstrun"-"$lastrun".root ]] ; then 
    echo "Runs_"$firstrun"_"$lastrun".root already exists!"
    exit
  fi 
  if [ -e $HIST_DIR/Sum_"$run".root ] ; then
  FILE_LIST=(${FILE_LIST[@]} "$HIST_DIR/Sum_"$run".root")
  fi
done

# sum all runs into single .root file
echo "::: Summing runs for $firstrun to $lastrun ..."
echo ""
echo "hadd -f $HIST_DIR/Runs_"$firstrun"-"$lastrun".root ${FILE_LIST[@]}"
  hadd -f $HIST_DIR/Runs_"$firstrun"-"$lastrun".root ${FILE_LIST[@]}
echo ""
echo "::: Summing runs for $firstrun to $lastrun ... [Done]"
