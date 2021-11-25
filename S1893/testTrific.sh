#!/bin/bash
analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis51994_000.root
#analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis51989_000.root
calFile=/tig/pterodon_data3/S1893/trificCal.cal
histDir=/tig/pterodon_data3/S1893/testHistograms
sortCode=/tig/pterodon_data3/S1893/SortCode/SortData


echo "Sorting TRIFIC data with the following code:"
echo "$sortCode $analyTree $calFile $histDir/trificHistogram.root"

$sortCode $analyTree $calFile $histDir/trificHistogram.root

