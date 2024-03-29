#!/bin/bash

declare -a FILE_LIST # array to hold file names
HIST_DIR=./HistFiles

firstrun=55585
lastrun=55588

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

# summing on a run by run basis
for run in `seq $firstrun $lastrun` ; do
  if [[ -e $HIST_DIR/Sum_$run.root ]] ; then 
    echo "Sum_$run.root already exists!"
    continue 
  fi 

  echo "::: Summing subruns for $run ... "
  echo ""
  for file in "$HIST_DIR"/Hist_"$run"_*
    do 
    FILE_LIST=(${FILE_LIST[@]} "$file")
    done

  # sum all subruns into single .root file
if [ ${#FILE_LIST[@]} -ne 0 ]; then
#echo "${#FILE_LIST[@]}   "
  echo "hadd -f $HIST_DIR/Sum_$run.root ${FILE_LIST[@]}"
  hadd -f $HIST_DIR/Sum_$run.root ${FILE_LIST[@]}

else
  if [ -e ${FILE_LIST[0]} ] ; then
  echo "cp ${FILE_LIST[0]} $HIST_DIR/Sum_$run.root"
  cp ${FILE_LIST[0]} $HIST_DIR/Sum_$run.root
  fi
fi
  # need to clear array after each loop 

  unset FILE_LIST
  echo ""
  echo "::: Summing subruns for $run ... [Done]"

done

# sum all individual runs together
echo "::: Summing All Runs :::"
  ./SumAllHist.sh $firstrun $lastrun
echo "::: [DONE] :::"
