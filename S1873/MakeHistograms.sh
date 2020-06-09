#!/bin/bash

# S1873 autosort script, originally by S. Gillespie

ANALDIR=AnalysisTrees
HISTDIR=HistFiles
DLEN=${#ANALDIR}
SORTCODE=SortCode/SortData
CFILE=SortCode/CalibrationFile.cal
CUTFILE=SortCode/cuts.root
TARGET=one # one, two, three
m=52046 #Only Sort AFTER this run 
n=52093 #Only Sort until this run 

if [ $# -eq 0 ]; then 
  echo "usage $0 <first run> <last run>"
  exit
fi

if [ $# -eq 1 ]; then 
  m=$1 - 1
  n=$1 + 1
fi

if [ $# -eq 2 ]; then
  m=$1 - 1
  n=$2 + 1
fi

if [ $# -eq 0 ]
	then

for f in $ANALDIR/*.root
 
do
# echo "$f"
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

else 

for f in $ANALDIR/*.root
 
do
# echo "$f"
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
