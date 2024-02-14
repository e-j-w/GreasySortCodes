CXXFLAGS = -std=c++0x -O2
GRSISORT = -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore

all: common.o SortDiagnosticsGRSI SortDiagnosticsSMOL DrawTIPWaveforms PrintRFData ComptonAsym_mca_SMOL ComptonAsym_mca_SMOL_EEGamma ComptonAsym_mca_SMOL_EEDopp ComptonAsym_mca_SMOL_EDopp ComptonAngle_SMOL ComptonAngle_SMOL_EDopp ComptonAngle_SMOL_EGammaEDopp RDCO_SMOL RDCO_SMOL_EDopp RDCO_SMOL_EGammaEDopp AngCorrSMOL AngDistSMOL SumSMOL CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EEGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EGamma_PIDsep_mca EDopp_PIDsep_mca TigEE_PIDsep_mca EGamma_mca_SMOL EGamma_noAB_mca_SMOL EGamma_no90deg_mca_SMOL EGamma_mca_SMOL_nopos15 EEGamma_mca_SMOL EEGamma_no90deg_mca_SMOL EEGamma_only180degcoinc_mca_SMOL EEGamma_only180degcoincnoAB_mca_SMOL EDopp_mca_SMOL EEDopp_mca_SMOL EEDopp_only180degcoinc_mca_SMOL EEDopp_mca_SMOL_ANDGate EEDopp_no90deg_mca_SMOL EDoppEGamma_mca_SMOL EDoppGateAngleEGamma_mca_SMOL EGammaEDopp_mca_SMOL ETIP_mca_SMOL TTGamma_SMOL TDoppTGamma_SMOL TTGamma_PIDsep TTipTGamma_PIDsep TIPRingPairs SeparatorTimingPID SeparatorSource

common.o: src/common.cxx src/position_arrays.cxx include/common.h include/evt_fmt.h
	g++ src/common.cxx -I./include $(CXXFLAGS) $(GRSISORT) -c -o common.o
SortDiagnosticsGRSI: src/SortDiagnosticsGRSI.cxx include/SortDiagnosticsGRSI.h common.o
	g++ src/SortDiagnosticsGRSI.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SortDiagnosticsGRSI
SortDiagnosticsSMOL: src/SortDiagnosticsSMOL.cxx include/SortDiagnosticsSMOL.h include/evt_fmt.h common.o
	g++ src/SortDiagnosticsSMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SortDiagnosticsSMOL
DrawTIPWaveforms: src/DrawTIPWaveforms.cxx
	g++ src/DrawTIPWaveforms.cxx -I./include $(CXXFLAGS) $(GRSISORT) -o DrawTIPWaveforms
PrintRFData: src/PrintRFData.cxx
	g++ src/PrintRFData.cxx -I./include $(CXXFLAGS) $(GRSISORT) -o PrintRFData
ComptonAsym_mca_SMOL: src/ComptonAsym_mca_SMOL.cxx include/ComptonAsym_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAsym_mca_SMOL
ComptonAsym_mca_SMOL_EEGamma: src/ComptonAsym_mca_SMOL_EEGamma.cxx include/ComptonAsym_mca_SMOL_EEGamma.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EEGamma.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAsym_mca_SMOL_EEGamma
ComptonAsym_mca_SMOL_EEDopp: src/ComptonAsym_mca_SMOL_EEDopp.cxx include/ComptonAsym_mca_SMOL_EEDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EEDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAsym_mca_SMOL_EEDopp
ComptonAsym_mca_SMOL_EDopp: src/ComptonAsym_mca_SMOL_EDopp.cxx include/ComptonAsym_mca_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAsym_mca_SMOL_EDopp
ComptonAngle_SMOL: src/ComptonAngle_SMOL.cxx include/ComptonAngle_SMOL.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAngle_SMOL
ComptonAngle_SMOL_EDopp: src/ComptonAngle_SMOL_EDopp.cxx include/ComptonAngle_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL_EDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAngle_SMOL_EDopp
ComptonAngle_SMOL_EGammaEDopp: src/ComptonAngle_SMOL_EGammaEDopp.cxx include/ComptonAngle_SMOL_EGammaEDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL_EGammaEDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ComptonAngle_SMOL_EGammaEDopp
RDCO_SMOL: src/RDCO_SMOL.cxx include/RDCO_SMOL.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o RDCO_SMOL
RDCO_SMOL_EDopp: src/RDCO_SMOL_EDopp.cxx include/RDCO_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL_EDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o RDCO_SMOL_EDopp
RDCO_SMOL_EGammaEDopp: src/RDCO_SMOL_EGammaEDopp.cxx include/RDCO_SMOL_EGammaEDopp.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL_EGammaEDopp.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o RDCO_SMOL_EGammaEDopp
AngCorrSMOL: src/AngCorrSMOL.cxx include/AngCorrSMOL.h include/evt_fmt.h common.o
	g++ src/AngCorrSMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o AngCorrSMOL
AngDistSMOL: src/AngDistSMOL.cxx include/AngDistSMOL.h include/evt_fmt.h common.o
	g++ src/AngDistSMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o AngDistSMOL
SumSMOL: src/SumSMOL.cxx include/SumSMOL.h include/evt_fmt.h common.o
	g++ src/SumSMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SumSMOL
CheckTimingWindows: src/CheckTimingWindows.cxx include/CheckTimingWindows.h common.o
	g++ src/CheckTimingWindows.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o CheckTimingWindows
PlotTimeSepPID: src/PlotTimeSepPID.cxx include/PlotTimeSepPID.h common.o
	g++ src/PlotTimeSepPID.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o PlotTimeSepPID
EGamma_PIDsep: src/EGamma_PIDsep.cxx include/EGamma_PIDsep.h common.o
	g++ src/EGamma_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_PIDsep
EEGamma_PIDsep: src/EEGamma_PIDsep.cxx include/EEGamma_PIDsep.h common.o
	g++ src/EEGamma_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEGamma_PIDsep
EDopp_PIDsep: src/EDopp_PIDsep.cxx include/EDopp_PIDsep.h common.o
	g++ src/EDopp_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_PIDsep
EEDopp_PIDsep: src/EEDopp_PIDsep.cxx include/EEDopp_PIDsep.h common.o
	g++ src/EEDopp_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_PIDsep
EDopp_ETIPtot_PIDsep: src/EDopp_ETIPtot_PIDsep.cxx include/EDopp_ETIPtot_PIDsep.h common.o
	g++ src/EDopp_ETIPtot_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_ETIPtot_PIDsep
EGamma_PIDsep_mca: src/EGamma_PIDsep_mca.cxx include/EGamma_PIDsep_mca.h common.o
	g++ src/EGamma_PIDsep_mca.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_PIDsep_mca
EDopp_PIDsep_mca: src/EDopp_PIDsep_mca.cxx include/EDopp_PIDsep_mca.h common.o
	g++ src/EDopp_PIDsep_mca.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_PIDsep_mca
TigEE_PIDsep_mca: src/TigEE_PIDsep_mca.cxx include/TigEE_PIDsep_mca.h common.o
	g++ src/TigEE_PIDsep_mca.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TigEE_PIDsep_mca
EGamma_mca_SMOL: src/EGamma_mca_SMOL.cxx include/EGamma_mca_SMOL.h common.o
	g++ src/EGamma_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_mca_SMOL
EGamma_noAB_mca_SMOL: src/EGamma_noAB_mca_SMOL.cxx include/EGamma_noAB_mca_SMOL.h common.o
	g++ src/EGamma_noAB_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_noAB_mca_SMOL
EGamma_no90deg_mca_SMOL: src/EGamma_no90deg_mca_SMOL.cxx include/EGamma_no90deg_mca_SMOL.h common.o
	g++ src/EGamma_no90deg_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_no90deg_mca_SMOL
EGamma_mca_SMOL_nopos15: src/EGamma_mca_SMOL_nopos15.cxx include/EGamma_mca_SMOL_nopos15.h common.o
	g++ src/EGamma_mca_SMOL_nopos15.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGamma_mca_SMOL_nopos15
EEGamma_mca_SMOL: src/EEGamma_mca_SMOL.cxx include/EEGamma_mca_SMOL.h common.o
	g++ src/EEGamma_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEGamma_mca_SMOL
EEGamma_no90deg_mca_SMOL: src/EEGamma_no90deg_mca_SMOL.cxx include/EEGamma_no90deg_mca_SMOL.h common.o
	g++ src/EEGamma_no90deg_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEGamma_no90deg_mca_SMOL
EEGamma_only180degcoinc_mca_SMOL: src/EEGamma_only180degcoinc_mca_SMOL.cxx include/EEGamma_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEGamma_only180degcoinc_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEGamma_only180degcoinc_mca_SMOL
EEGamma_only180degcoincnoAB_mca_SMOL: src/EEGamma_only180degcoincnoAB_mca_SMOL.cxx include/EEGamma_only180degcoincnoAB_mca_SMOL.h common.o
	g++ src/EEGamma_only180degcoincnoAB_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEGamma_only180degcoincnoAB_mca_SMOL
EDopp_mca_SMOL: src/EDopp_mca_SMOL.cxx include/EDopp_mca_SMOL.h common.o
	g++ src/EDopp_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDopp_mca_SMOL
EEDopp_mca_SMOL: src/EEDopp_mca_SMOL.cxx include/EEDopp_mca_SMOL.h common.o
	g++ src/EEDopp_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_mca_SMOL
EEDopp_only180degcoinc_mca_SMOL: src/EEDopp_only180degcoinc_mca_SMOL.cxx include/EEDopp_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEDopp_only180degcoinc_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_only180degcoinc_mca_SMOL
EEDopp_mca_SMOL_ANDGate: src/EEDopp_mca_SMOL_ANDGate.cxx include/EEDopp_mca_SMOL_ANDGate.h common.o
	g++ src/EEDopp_mca_SMOL_ANDGate.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_mca_SMOL_ANDGate
EEDopp_no90deg_mca_SMOL: src/EEDopp_no90deg_mca_SMOL.cxx include/EEDopp_no90deg_mca_SMOL.h common.o
	g++ src/EEDopp_no90deg_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EEDopp_no90deg_mca_SMOL
EDoppEGamma_mca_SMOL: src/EDoppEGamma_mca_SMOL.cxx include/EDoppEGamma_mca_SMOL.h common.o
	g++ src/EDoppEGamma_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDoppEGamma_mca_SMOL
EDoppGateAngleEGamma_mca_SMOL: src/EDoppGateAngleEGamma_mca_SMOL.cxx include/EDoppGateAngleEGamma_mca_SMOL.h common.o
	g++ src/EDoppGateAngleEGamma_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EDoppGateAngleEGamma_mca_SMOL
EGammaEDopp_mca_SMOL: src/EGammaEDopp_mca_SMOL.cxx include/EGammaEDopp_mca_SMOL.h common.o
	g++ src/EGammaEDopp_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o EGammaEDopp_mca_SMOL
ETIP_mca_SMOL: src/ETIP_mca_SMOL.cxx include/ETIP_mca_SMOL.h common.o
	g++ src/ETIP_mca_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o ETIP_mca_SMOL
TTGamma_SMOL: src/TTGamma_SMOL.cxx include/TTGamma_SMOL.h common.o
	g++ src/TTGamma_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TTGamma_SMOL
TDoppTGamma_SMOL: src/TDoppTGamma_SMOL.cxx include/TDoppTGamma_SMOL.h common.o
	g++ src/TDoppTGamma_SMOL.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TDoppTGamma_SMOL
TTGamma_PIDsep: src/TTGamma_PIDsep.cxx include/TTGamma_PIDsep.h common.o
	g++ src/TTGamma_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TTGamma_PIDsep
TTipTGamma_PIDsep: src/TTipTGamma_PIDsep.cxx include/TTipTGamma_PIDsep.h common.o
	g++ src/TTipTGamma_PIDsep.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TTipTGamma_PIDsep
TIPRingPairs: src/TIPRingPairs.cxx include/TIPRingPairs.h common.o
	g++ src/TIPRingPairs.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o TIPRingPairs
SeparatorTimingPID: src/SeparatorTimingPID.cxx include/SeparatorTimingPID.h include/evt_fmt.h common.o
	g++ src/SeparatorTimingPID.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SeparatorTimingPID
SeparatorSource: src/SeparatorSource.cxx include/SeparatorSource.h include/evt_fmt.h common.o
	g++ src/SeparatorSource.cxx -I./include $(CXXFLAGS) $(GRSISORT) common.o -o SeparatorSource
clean:
	rm -rf *~ SortDiagnosticsGRSI SortDiagnosticsSMOL DrawTIPWaveforms PrintRFData ComptonAsym_mca_SMOL ComptonAsym_mca_SMOL_EEGamma ComptonAsym_mca_SMOL_EEDopp ComptonAsym_mca_SMOL_EDopp ComptonAngle_SMOL ComptonAngle_SMOL_EDopp ComptonAngle_SMOL_EGammaEDopp RDCO_SMOL RDCO_SMOL_EDopp RDCO_SMOL_EGammaEDopp AngCorrSMOL AngDistSMOL SumSMOL CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EEGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EGamma_PIDsep_mca EDopp_PIDsep_mca TigEE_PIDsep_mca EGamma_mca_SMOL EGamma_noAB_mca_SMOL EGamma_no90deg_mca_SMOL EEGamma_mca_SMOL EEGamma_no90deg_mca_SMOL EEGamma_only180degcoinc_mca_SMOL EEGamma_only180degcoincnoAB_mca_SMOL EDopp_mca_SMOL EGamma_mca_SMOL_nopos15 EEDopp_mca_SMOL EEDopp_only180degcoinc_mca_SMOL EEDopp_mca_SMOL_ANDGate EEDopp_no90deg_mca_SMOL EDoppEGamma_mca_SMOL EDoppGateAngleEGamma_mca_SMOL EGammaEDopp_mca_SMOL ETIP_mca_SMOL TTGamma_SMOL TDoppTGamma_SMOL TTGamma_PIDsep TTipTGamma_PIDsep TIPRingPairs SeparatorTimingPID SeparatorSource common.o
