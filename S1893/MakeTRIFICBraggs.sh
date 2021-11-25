#!/bin/bash

grsisort -l

auto tree = "/tig/pterodon_data3/S1893/testAnalysisTrees/analysis$@_000.root"

.x SortCode/TrificTestSortv2.C
