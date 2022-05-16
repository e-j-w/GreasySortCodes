CXXFLAGS = -std=c++0x -O2
GRSISORT = -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lMinuit -lPyROOT

all: common.o SortCode DrawTIPWaveforms CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EDopp_PIDsep_mca

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
EDopp_PIDsep: src/EDopp_PIDsep.cxx include/EDopp_PIDsep.h common.o
	g++ src/EDopp_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_PIDsep
EEDopp_PIDsep: src/EEDopp_PIDsep.cxx include/EEDopp_PIDsep.h common.o
	g++ src/EEDopp_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_PIDsep
EDopp_ETIPtot_PIDsep: src/EDopp_ETIPtot_PIDsep.cxx include/EDopp_ETIPtot_PIDsep.h common.o
	g++ src/EDopp_ETIPtot_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_ETIPtot_PIDsep
EDopp_PIDsep_mca: src/EDopp_PIDsep_mca.cxx include/EDopp_PIDsep_mca.h common.o
	g++ src/EDopp_PIDsep_mca.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_PIDsep_mca
clean:
	rm -rf *~ SortCode DrawTIPWaveforms CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EDopp_PIDsep_mca common.o
