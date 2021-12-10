#!/bin/bash

ADIR=AnalysisTrees
HISTDIR=HistFiles
DLEN=${#ADIR}
SORTCODE=SortCode/SortData
CFILE=CalibrationFile.cal
#CFILE=CalibrationFileNoCal.cal
TARGET=one
m=55777 #Only Sort after this run
n=80000 #Only Sort until this run 

if [ ! -d $HISTDIR ]; then
 mkdir $HISTDIR
fi


if [ $# -eq 0 ]
	then

for f in $ADIR/*.root
 
do
# echo "$f"[S1801@pterodon S1801]$ ./MakeHis
 g=${f:DLEN+9}
 h=${g:0:${#g}-5} 
 i=${g:0:${#g}-9}
 echo "$i"

 HFILE=$HISTDIR/Hist_$h.root

 if [ $i -gt $m ] && [ $n -gt $i ];
  then
   echo "$SORTCODE $f $CFILE $HFILE"
   $SORTCODE $f $CFILE $HFILE
  fi
done

else 

for f in $ADIR/analysis"$@"_*.root;
 
do
 echo "$f"
 g=${f:DLEN+9}
 h=${g:0:${#g}-5} 
 i=${g:0:${#g}-9} 

 HFILE=$HISTDIR/Hist_$h.root

 if [ $i -gt $m ] && [ $n -gt $i ];
  then
   echo "$SORTCODE $f $CFILE $HFILE $TARGET"
   $SORTCODE $f $CFILE $HFILE $TARGET
  fi
done

fi
