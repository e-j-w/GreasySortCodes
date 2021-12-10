#!/bin/bash

ANALDIR=../AnalysisTrees
HISTDIR=../HistFiles/Hist_Matt
DLEN=${#ANALDIR}
SORTCODE=./SortData_Matt
CFILE=../CalibrationFile.cal
CUTFILE=SortCode/cuts.root
TARGET=one # one, two, three
m=55571 #Only Sort after this run 
n=99999 #Only Sort until this run
if [ $# -eq 0 ]
	then

for f in $ANALDIR/*.root
 
do
# echo "$f"
 g=${f:DLEN+9}
 h=${g:0:${#g}-5} 
 i=${g:0:${#g}-9} 

 HFILE=$HISTDIR/Hist_Matt_$h.root

 if [ $i -gt $m ] && [ $n -gt $i ];
  then
  if [ ! -f $HFILE ];
   then
    if [ ! -f $CUTFILE ];
     then
      echo "$SORTCODE $f $CFILE $HFILE $TARGET"
      $SORTCODE $f $CFILE $HFILE $TARGET
     else
      echo "$SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE"
      $SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE    
     fi 
   fi
   if [ -f $HFILE ]
    then
     if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $f ];
      then
       if [ ! -f $CUTFILE ];
        then
         echo "$SORTCODE $f $CFILE $HFILE $TARGET"
         $SORTCODE $f $CFILE $HFILE $TARGET
        else
         echo "$SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE"
         $SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE    
       fi      
     fi
    fi
  fi
done

else 

for f in $ANALDIR/analysis"$@"_*.root;
 
do
 echo "$f"
 g=${f:DLEN+9}
 h=${g:0:${#g}-5} 
 i=${g:0:${#g}-9} 

 HFILE=$HISTDIR/Hist_$h.root

 if [ $i -gt $m ] && [ $n -gt $i ];
  then
  if [ ! -f $HFILE ];
   then
    if [ ! -f $CUTFILE ];
     then
      echo "$SORTCODE $f $CFILE $HFILE $TARGET"
      $SORTCODE $f $CFILE $HFILE $TARGET
     else
      echo "$SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE"
      $SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE    
     fi 
   fi
   if [ -f $HFILE ]
    then
     if [ $HFILE -ot $SORTCODE ] || [ $HFILE -ot $f ];
      then
       if [ ! -f $CUTFILE ];
        then
         echo "$SORTCODE $f $CFILE $HFILE $TARGET"
         $SORTCODE $f $CFILE $HFILE $TARGET
        else
         echo "$SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE"
         $SORTCODE $f $CFILE $HFILE $TARGET $CUTFILE    
       fi      
     fi
    fi
  fi
done

fi
