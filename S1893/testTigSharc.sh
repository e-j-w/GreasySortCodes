#!/bin/bash
#analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis51994_000.root
#analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis29035_000.root
analyTree=/tig/pterodon_data3/S1893/testAnalysisTrees/analysis55380_000.root
#calFile=/tig/pterodon_data3/S1893/S1389Cal.cal
calFile=/tig/pterodon_data3/S1893/CalibrationFile.cal
histDir=/tig/pterodon_data3/S1893/testHistograms
sortCode=/tig/pterodon_data3/S1893/SortCode/SortData


echo "Sorting Tigress-Sharc data test alpha with the following code:"
echo "$sortCode $analyTree $calFile $histDir/tigSharcHistogram.root"

$sortCode $analyTree $calFile $histDir/tigSharcHistogram.root

