#!/bin/bash

ANALDIR=AnalysisTrees
DLEN=${#ANALDIR}
SORTCODE=SortCode/Clock
#CFILE=SortCode/CalibrationFile.cal
CFILE=CalibrationFile.cal

for f in $ANALDIR/analysis"$@"_000.root
 
do
 echo "$f"
 g=${f:DLEN+9}
 h=${g:0:${#g}-5} 
 HFILE=Clock.root

 echo "$SORTCODE $f $CFILE $HFILE"
 $SORTCODE $f $CFILE $HFILE
done

