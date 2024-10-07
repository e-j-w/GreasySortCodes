CXXFLAGS = -O2 -Wl,--copy-dt-needed-entries
GRSISORT = $(shell grsi-config --cflags --all-libs --GRSIData-libs) -I$(GRSISYS)/GRSIData/include $(shell root-config --cflags --libs)
#GRSISORT = -I/data1/jwilliams/Code/GRSISort-centos7/include -L/data1/jwilliams/Code/GRSISort-centos7/lib -lTDetector -lTGRSIint -lTFormat -lTGRSIint -lTFormat -lTLoops -lTKinematics -lTSRIM -lTBetaDecay -lTCal -lTReaction -lTPeakFitting -lTGRSIFrame -lTSuppressed -lTBgo -lTPulseAnalyzer -lGROOT -lTGRSIFit -lTNucleus -lTRawFile -lTDataParser -lTHistogramming -lTGUI -lTLevelScheme -lXMLParser -lXMLIO -lX11 -lXpm -lProof -lGuiHtml -lMinuit -lSpectrum -lMathMore -L/data1/jwilliams/Code/GRSISort-centos7/GRSIData/lib -lTAngularCorrelation -lTCSM -lTDescant -lTGRSIDataParser -lTGRSIFormat -lTGenericDetector -lTGriffin -lTLaBr -lTMidas -lTPaces -lTRF -lTS3 -lTSceptar -lTSharc -lTSiLi -lTTAC -lTTigress -lTTriFoil -lTZeroDegree -lTEmma -lTTrific -I/data1/jwilliams/Code/GRSISort-centos7/GRSIData/include -pthread -std=c++17 -m64 -I/data1/jwilliams/Code/root_6.32.02/include -L/data1/jwilliams/Code/root_6.32.02/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/data1/jwilliams/Code/root_6.32.02/lib -pthread -lm -ldl -rdynamic -lTreePlayer -lMathMore -lX11

all: common.o SortDiagnosticsSMOL SeparatorSource

common.o: src/common.cxx src/position_arrays.cxx include/common.h include/evt_fmt.h
	g++ src/common.cxx -I./include $(CXXFLAGS) $(GRSISORT) -c -o common.o
SortDiagnosticsSMOL: src/SortDiagnosticsSMOL.cxx include/SortDiagnosticsSMOL.h include/evt_fmt.h common.o
	g++ src/SortDiagnosticsSMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) -o SortDiagnosticsSMOL
SeparatorSource: src/SeparatorSource.cxx include/SeparatorSource.h include/evt_fmt.h common.o
	g++ src/SeparatorSource.cxx -I./include common.o $(CXXFLAGS) $(GRSISORT) -o SeparatorSource
clean:
	rm -rf *~ SortDiagnosticsSMOL SeparatorSource common.o
