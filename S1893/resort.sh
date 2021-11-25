#!/bin/bash

for i in 55713; do
  echo $i
  /bin/bash ./MakeAnalysis.sh $i
  /bin/bash ./MakeHistograms.sh $i
done
