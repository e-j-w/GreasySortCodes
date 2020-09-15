#!/bin/bash

# GRSIsort analysis tree helper script, original by S. Gillespie

DATADIR=MidasFiles
FRAGDIR=FragmentTrees
ANALDIR=AnalysisTrees

DLEN=${#DATADIR}
 
m=51561

if [ $# -eq 0 ]
	then
	for f in $DATADIR/*.mid; do
	 g=${f:DLEN+4}
	 h=${g:0:${#g}-4} 
	 i=${g:0:${#g}-8} 
	 FFILE=$FRAGDIR/fragment$h.root
	 AFILE=$ANALDIR/analysis$h.root
	 CFILE=SortCode/CalibrationFile.cal
	 #echo "Processing run$g $h $i"

		if [ ! -f $AFILE ] && [ $i -gt $m ];  
				then
					echo "File $AFILE does not exist."
					grsisort -laq $f $CFILE --suppress-errors --write-frag-tree --word-count-offset=-1 --sort-depth 1000000000 --build-window 5000
					mv -f analysis$h.root $ANALDIR
					mv -f fragment$h.root $FRAGDIR

		fi
			 
		if [ -f $AFILE ]
		then 

			if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
					then
						echo "File $AFILE exists but is older than $f"
						grsisort -laq $f $CFILE --suppress-errors --write-frag-tree --word-count-offset=-1 --sort-depth 1000000000 --build-window 5000
						mv -f analysis$h.root $ANALDIR
						mv -f fragment$h.root $FRAGDIR
			fi
		fi
	done 

else 
	for f in $DATADIR/run"$@"_*.mid;

	do
	 g=${f:DLEN+4}
	 h=${g:0:${#g}-4} 
	 i=${g:0:${#g}-8} 
	 FFILE=$FRAGDIR/fragment$h.root
	 AFILE=$ANALDIR/analysis$h.root
	 CFILE=SortCode/CalibrationFile.cal
	 echo "Processing run$g "

		if [ ! -f $AFILE ];  
				then
					echo "File $AFILE does not exist."
					grsisort -laq $f $CFILE --suppress-errors --write-frag-tree --word-count-offset=-1 --sort-depth 1000000000 --build-window 5000
					mv -f analysis$h.root $ANALDIR
					mv -f fragment$h.root $FRAGDIR

		fi

		if [ -f $AFILE ]
		then 
			
			if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
					then
						echo "File $AFILE exists but is older than $f"
						grsisort -laq $f $CFILE --suppress-errors --write-frag-tree --word-count-offset=-1 --sort-depth 1000000000 --build-window 5000
						mv -f analysis$h.root $ANALDIR
						mv -f fragment$h.root $FRAGDIR
			fi
		fi
	done
fi
