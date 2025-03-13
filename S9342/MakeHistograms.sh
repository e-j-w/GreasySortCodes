#!/bin/bash

ADIR=AnalysisTrees
SMOLDIR=SmolFiles
HISTDIR=HistFiles
DLEN=${#ADIR}
SORTCODE=SortDiagnosticsSMOL
CFILE=GRIFFIN-Cal-File-Run29604.cal
FIRSTSUBRUN=0
LASTSUBRUN=10000000

if [ ! -d $HISTDIR ]; then
mkdir $HISTDIR
fi
if [ ! -d $SMOLDIR ]; then
mkdir $SMOLDIR
fi

if [ "$#" -lt 1 ]
then
	echo "Provide a run number as an argument!"
	echo "First and last subruns can be specified as 2nd and 3rd arguments (optional)."
else

	if [ $1 -ge 29711 ];
	then
		CFILE=GRIFFIN-Cal-File-Run29711.cal
	elif [ $1 -ge 29670 ];
	then
		CFILE=GRIFFIN-Cal-File-Run29671.cal
	fi
	
	SMOLFILE=$SMOLDIR/run"$1".smole6
	HFILE=$HISTDIR/Hist_"$1".root
	if [ "$#" -eq 3 ]
	then
		FIRSTSUBRUN=$2
		LASTSUBRUN=$3
		SMOLFILE=$SMOLDIR/run"$1"_"$FIRSTSUBRUN"_"$LASTSUBRUN".smole6
		HFILE=$HISTDIR/Hist_"$1"_"$FIRSTSUBRUN"_"$LASTSUBRUN".root
	fi
	if [ ! -f $SMOLFILE ];
	then
		# empty the list file
		> run"$1".list
		for f in $ADIR/analysis"$1"*.root
		do
			g=${f:DLEN+9}
			h=${g:0:${#g}-5} 
			i=${g:6:${#g}-11}
			#echo "$i $FIRSTSUBRUN $LASTSUBRUN"

			if [ $i -ge $FIRSTSUBRUN ] && [ $i -le $LASTSUBRUN ];
			then
				echo "$f"
				ls $f >> run"$1".list
			fi

		done
		./E6Sort/SeparatorSource run"$1".list $CFILE $SMOLFILE
		rm run"$1".list
	fi
	./E6Sort/$SORTCODE $SMOLFILE $HFILE
fi
