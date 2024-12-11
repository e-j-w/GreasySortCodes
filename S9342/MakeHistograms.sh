#!/bin/bash

ADIR=AnalysisTrees
SMOLDIR=SmolFiles
HISTDIR=HistFiles
DLEN=${#ADIR}
SORTCODE=SortDiagnosticsSMOL
CFILE=CalibrationFileGRIFFIN.cal
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
	SMOLFILE=$SMOLDIR/run"$1".smole6
	HFILE=$HISTDIR/Hist_"$1".root
	if [ "$#" -eq 3 ]
	then
		FIRSTSUBRUN=$2
		LASTSUBRUN=$3
		SMOLFILE=$SMOLDIR/run"$1"_"$FIRSTSUBRUN"_"$LASTSUBRUN".smole6
		HFILE=$HISTDIR/Hist_"$1"_"$FIRSTSUBRUN"_"$LASTSUBRUN".root
	fi
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
	if [ ! -f $SMOLFILE ];
	then
		./E6Sort/SeparatorSource run"$1".list $CFILE $SMOLFILE
	fi
	./E6Sort/$SORTCODE $SMOLFILE $HFILE
	rm run"$1".list
fi
