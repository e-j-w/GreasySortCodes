#!/bin/bash

# GRSIsort analysis tree helper script, original by S. Gillespie with modifications by J. Williams

#DATADIR=MidasFiles
DATADIR=/tig/grifstore1/grifalt/schedule147/S9342
#DATADIR=midas-92mNb
ADIR=AnalysisTrees
GRSIARGS="--suppress-errors --recommended --sort-depth=10000000"

if [ ! -d $ADIR ]; then
	mkdir $ADIR
fi

DLEN=${#DATADIR}
 
m=10000
CFILE=GRIFFIN-Cal-File-Run29604.cal

if [ $# -eq 1 ]
	then
	
	if [ $1 -ge 29711 ]
	then
		CFILE=GRIFFIN-Cal-File-Run29711.cal
	elif [ $1 -ge 29670 ]
	then
		CFILE=GRIFFIN-Cal-File-Run29671.cal
	fi
	
	for f in $DATADIR/run"$@"_*.mid;

	do
	 g=${f:DLEN+4}
	 h=${g:0:${#g}-4} 
	 i=${g:0:${#g}-8}
	 AFILE=$ADIR/analysis$h.root
	 echo "Processing run$g "

		if [ ! -f $AFILE ];  
		then
			echo "File $AFILE does not exist."
			grsisort -laq $f $CFILE $GRSIARGS
			mv -f analysis$h.root $ADIR
			rm -f fragment$h.root

		fi

		if [ -f $AFILE ]
		then
			if [ "$AFILE" -ot "$f" ] && [ $i -gt $m ];
					then
						echo "File $AFILE exists but is older than $f"
						grsisort -laq $f $CFILE $GRSIARGS
						mv -f analysis$h.root $ADIR
						rm -f fragment$h.root
			fi
		fi
	done

else 
	echo "Give a run number as an argument."
fi
