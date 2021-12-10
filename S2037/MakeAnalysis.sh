#!/bin/bash
#DATADIR=/tig/alphadon_data1/LabTests
DATADIR=/tig/tigstore01/Schedule141/S2037
ADIR=AnalysisTrees
FDIR=FragmentTrees
CFILE=CalibrationFile.cal
#CFILE=CalibrationFileClean.cal
DLEN=${#DATADIR}
y=$(printf "%05d" $@)

if [ ! -d $FDIR ]; then
 mkdir $FDIR
fi
if [ ! -d $ADIR ]; then
 mkdir $ADIR
fi

for f in $DATADIR/run"$y"_*.mid; # Sorts sub runs from a single file

do
 g=${f:DLEN+4}
 h=${g:0:${#g}-4} 
 i=${g:0:${#g}-8} 
 FFILE=fragment$h.root
 AFILE=analysis$h.root
 echo "Processing run$g "

 if [ ! -f $AFILE ];  
  then
  echo "File $AFILE does not exist."
    grsisort --recommended $f $CFILE --word-count-offset=0 --sort-depth 100000000
    #grsisort --recommended $f $CFILE --build-window 100000 --word-count-offset=0 --sort-depth 1000000000

 else 
   if [ -f $AFILE ]
   then 
    if [ "$AFILE" -ot "$f" ];
     then
     echo "File $AFILE exists but is older than $f"
      grsisort --recommended $f $CFILE --word-count-offset=0 --sort-depth 1000000000
      #grsisort --recommended $f $CFILE --build-window 100000 --word-count-offset=0 --sort-depth 1000000000
    fi
   fi
 fi

mv $FFILE $FDIR
mv $AFILE $ADIR

done





