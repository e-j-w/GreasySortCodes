#!/bin/bash

# GRSIsort analysis tree helper script, original by S. Gillespie

#DATADIR=MidasFiles
DATADIR=/tig/grifstore1/grifalt/schedule147/S2231_S2196_Nov2024
FRAGDIR=FragmentTrees
ADIR=AnalysisTrees

if [ ! -d $ADIR ]; then
	mkdir $ADIR
fi

DLEN=${#DATADIR}
 
m=10000
CFILE=CalibrationFileGRIFFIN.cal

if [ $# -eq 0 ]
	then
	for f in $DATADIR/*.mid; do
	 g=${f:DLEN+4}
	 h=${g:0:${#g}-4} 
	 i=${g:0:${#g}-8} 
	 FFILE=$FRAGDIR/fragment$h.root
	 AFILE=$ADIR/analysis$h.root
	 
	 #echo "Processing run$g $h $i"

		if [ ! -f $AFILE ] && [ $i -gt $m ];  
				then
					echo "File $AFILE does not exist."
					grsisort -laq $f $CFILE --suppress-errors --recommended --sort-depth=10000000
					mv -f analysis$h.root $ADIR
					mv -f fragment$h.root $FRAGDIR

		fi
			 
		if [ -f $AFILE ]
		then 

			if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
					then
						echo "File $AFILE exists but is older than $f"
						grsisort -laq $f $CFILE --suppress-errors --recommended --sort-depth=10000000
						mv -f analysis$h.root $ADIR
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
	 AFILE=$ADIR/analysis$h.root
	 echo "Processing run$g "

		if [ ! -f $AFILE ];  
				then
					echo "File $AFILE does not exist."
					grsisort -laq $f $CFILE --suppress-errors --recommended --sort-depth=10000000
					mv -f analysis$h.root $ADIR
					mv -f fragment$h.root $FRAGDIR

		fi

		if [ -f $AFILE ]
		then
			if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
					then
						echo "File $AFILE exists but is older than $f"
						grsisort -laq $f $CFILE --suppress-errors --recommended --sort-depth=10000000
						mv -f analysis$h.root $ADIR
						mv -f fragment$h.root $FRAGDIR
			fi
		fi
	done
fi
