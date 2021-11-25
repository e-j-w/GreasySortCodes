#!/bin/bash
analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis54545_000.root
calFile=/tig/pterodon_data3/S1893/sharcCal.cal
histDir=/tig/pterodon_data3/S1893/testHistograms
sortCode=/tig/pterodon_data3/S1893/SortCode/SortData

echo "Sorting SHARC data with the following code:"
echo "$sortCode $analyTree $calFile $histDir/sharcHistogram.root"

$sortCode $analyTree $calFile $histDir/sharcHistogram.root

