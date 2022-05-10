CXXFLAGS = -std=c++0x -O2
GRSISORT = -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT

all: common.o SortCode DrawTIPWaveforms CheckTimingWindows PlotTimeSepPID EGamma_PIDsep

common.o: src/common.cxx include/common.h
	g++ src/common.cxx -I./include $(CXXFLAGS) $(GRSISORT) -c -o common.o
SortCode: src/SortCode.cxx include/SortCode.h common.o
	g++ src/SortCode.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SortCode 
DrawTIPWaveforms: src/DrawTIPWaveforms.cxx
	g++ src/DrawTIPWaveforms.cxx -I./include $(CXXFLAGS) $(GRSISORT) -o DrawTIPWaveforms
CheckTimingWindows: src/CheckTimingWindows.cxx include/CheckTimingWindows.h common.o
	g++ src/CheckTimingWindows.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o CheckTimingWindows
PlotTimeSepPID: src/PlotTimeSepPID.cxx include/PlotTimeSepPID.h common.o
	g++ src/PlotTimeSepPID.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o PlotTimeSepPID
EGamma_PIDsep: src/EGamma_PIDsep.cxx include/EGamma_PIDsep.h common.o
	g++ src/EGamma_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_PIDsep
clean:
	rm -rf *~ SortCode DrawTIPWaveforms CheckTimingWindows PlotTimeSepPID EGamma_PIDsep common.o
