LINKFLAGS = -Wl,--no-as-needed $(shell root-config --cflags --libs --glibs) -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm
CXXFLAGS = -O2 -Wl,--copy-dt-needed-entries
GRSISORT = -L/opt/local/lib -lX11 -lXpm $(shell grsi-config --cflags --all-libs --GRSIData-libs) -I$(GRSISYS)/GRSIData/include
#GRSISORT = -I/data1/jwilliams/Code/GRSISort-centos7/include -L/data1/jwilliams/Code/GRSISort-centos7/lib -lTDetector -lTGRSIint -lTFormat -lTGRSIint -lTFormat -lTLoops -lTKinematics -lTSRIM -lTBetaDecay -lTCal -lTReaction -lTPeakFitting -lTGRSIFrame -lTSuppressed -lTBgo -lTPulseAnalyzer -lGROOT -lTGRSIFit -lTNucleus -lTRawFile -lTDataParser -lTHistogramming -lTGUI -lTLevelScheme -lXMLParser -lXMLIO -lX11 -lXpm -lProof -lGuiHtml -lMinuit -lSpectrum -lMathMore -L/data1/jwilliams/Code/GRSISort-centos7/GRSIData/lib -lTAngularCorrelation -lTCSM -lTDescant -lTGRSIDataParser -lTGRSIFormat -lTGenericDetector -lTGriffin -lTLaBr -lTMidas -lTPaces -lTRF -lTS3 -lTSceptar -lTSharc -lTSiLi -lTTAC -lTTigress -lTTip -lTTriFoil -lTZeroDegree -lTEmma -lTTrific -I/data1/jwilliams/Code/GRSISort-centos7/GRSIData/include -pthread -std=c++17 -m64 -I/data1/jwilliams/Code/root_6.32.02/include -L/data1/jwilliams/Code/root_6.32.02/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/data1/jwilliams/Code/root_6.32.02/lib -pthread -lm -ldl -rdynamic -lTreePlayer -lMathMore -lX11

all: common.o SortDiagnosticsGRSI SortDiagnosticsSMOL DrawTIPWaveforms PrintRFData ComptonAsym_mca_SMOL ComptonAsym_mca_SMOL_EEGamma ComptonAsym_mca_SMOL_EEDopp ComptonAsym_mca_SMOL_EDopp ComptonAngle_SMOL ComptonAngle_SMOL_EDopp ComptonAngle_SMOL_EGammaEDopp RDCO_SMOL RDCO_SMOL_EDopp RDCO_SMOL_EGammaEDopp AngCorrSMOL AngDistSMOL SumSMOL CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EEGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EGamma_PIDsep_mca EDopp_PIDsep_mca TigEE_PIDsep_mca EGamma_mca_SMOL EGamma_noAB_mca_SMOL EGamma_no90deg_mca_SMOL EGamma_mca_SMOL_nopos15 EGamma_phi_mca_SMOL EEGamma_mca_SMOL EEGamma_mult2timing_mca_SMOL EEGamma_noAB_mult2timing_mca_SMOL EEGamma_noAB_noSupp_mult2timing_mca_SMOL EEGamma_ANDgate_mca_SMOL EEGamma_no90deg_mca_SMOL EEGamma_only180degcoinc_mca_SMOL EEGamma_only180degcoincnoAB_mca_SMOL EDopp_mca_SMOL EDopp_noAB_mca_SMOL EDopp_noAB_mult2up_mca_SMOL EDopp_phi_mca_SMOL EEDopp_mca_SMOL EEDopp_noAB_mca_SMOL EGamma_only180degcoinc_mca_SMOL EDopp_only180degcoinc_mca_SMOL EEDopp_only180degcoinc_mca_SMOL EEEGamma_only180degcoinc_mca_SMOL EEEDopp_only180degcoinc_mca_SMOL EEDopp_ANDGate_mca_SMOL EEDopp_ORGate_mca_SMOL EEDopp_ORGate_noAB_mca_SMOL EEDopp_no90deg_mca_SMOL EDoppEGamma_mca_SMOL EDoppGateAngleEGamma_mca_SMOL EGammaEDopp_mca_SMOL EGammaEDopp_noAB_mca_SMOL ETIP_mca_SMOL TTGamma_SMOL TDoppTGamma_SMOL TTGamma_PIDsep TTipTGamma_PIDsep TIPRingPairs SeparatorTimingPID SeparatorSource

common.o: src/common.cxx src/position_arrays.cxx include/common.h include/evt_fmt.h
	g++ src/common.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -c -o common.o
SortDiagnosticsGRSI: src/SortDiagnosticsGRSI.cxx include/SortDiagnosticsGRSI.h common.o
	g++ src/SortDiagnosticsGRSI.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SortDiagnosticsGRSI
SortDiagnosticsSMOL: src/SortDiagnosticsSMOL.cxx include/SortDiagnosticsSMOL.h include/evt_fmt.h common.o
	g++ src/SortDiagnosticsSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SortDiagnosticsSMOL
DrawTIPWaveforms: src/DrawTIPWaveforms.cxx
	g++ src/DrawTIPWaveforms.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o DrawTIPWaveforms
PrintRFData: src/PrintRFData.cxx
	g++ src/PrintRFData.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o PrintRFData
ComptonAsym_mca_SMOL: src/ComptonAsym_mca_SMOL.cxx include/ComptonAsym_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAsym_mca_SMOL
ComptonAsym_mca_SMOL_EEGamma: src/ComptonAsym_mca_SMOL_EEGamma.cxx include/ComptonAsym_mca_SMOL_EEGamma.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EEGamma.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAsym_mca_SMOL_EEGamma
ComptonAsym_mca_SMOL_EEDopp: src/ComptonAsym_mca_SMOL_EEDopp.cxx include/ComptonAsym_mca_SMOL_EEDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EEDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAsym_mca_SMOL_EEDopp
ComptonAsym_mca_SMOL_EDopp: src/ComptonAsym_mca_SMOL_EDopp.cxx include/ComptonAsym_mca_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAsym_mca_SMOL_EDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAsym_mca_SMOL_EDopp
ComptonAngle_SMOL: src/ComptonAngle_SMOL.cxx include/ComptonAngle_SMOL.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAngle_SMOL
ComptonAngle_SMOL_EDopp: src/ComptonAngle_SMOL_EDopp.cxx include/ComptonAngle_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL_EDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAngle_SMOL_EDopp
ComptonAngle_SMOL_EGammaEDopp: src/ComptonAngle_SMOL_EGammaEDopp.cxx include/ComptonAngle_SMOL_EGammaEDopp.h include/evt_fmt.h common.o
	g++ src/ComptonAngle_SMOL_EGammaEDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ComptonAngle_SMOL_EGammaEDopp
RDCO_SMOL: src/RDCO_SMOL.cxx include/RDCO_SMOL.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o RDCO_SMOL
RDCO_SMOL_EDopp: src/RDCO_SMOL_EDopp.cxx include/RDCO_SMOL_EDopp.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL_EDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o RDCO_SMOL_EDopp
RDCO_SMOL_EGammaEDopp: src/RDCO_SMOL_EGammaEDopp.cxx include/RDCO_SMOL_EGammaEDopp.h include/evt_fmt.h common.o
	g++ src/RDCO_SMOL_EGammaEDopp.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o RDCO_SMOL_EGammaEDopp
AngCorrSMOL: src/AngCorrSMOL.cxx include/AngCorrSMOL.h include/evt_fmt.h common.o
	g++ src/AngCorrSMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o AngCorrSMOL
AngDistSMOL: src/AngDistSMOL.cxx include/AngDistSMOL.h include/evt_fmt.h common.o
	g++ src/AngDistSMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o AngDistSMOL
SumSMOL: src/SumSMOL.cxx include/SumSMOL.h include/evt_fmt.h common.o
	g++ src/SumSMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SumSMOL
CheckTimingWindows: src/CheckTimingWindows.cxx include/CheckTimingWindows.h common.o
	g++ src/CheckTimingWindows.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o CheckTimingWindows
PlotTimeSepPID: src/PlotTimeSepPID.cxx include/PlotTimeSepPID.h common.o
	g++ src/PlotTimeSepPID.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o PlotTimeSepPID
EGamma_PIDsep: src/EGamma_PIDsep.cxx include/EGamma_PIDsep.h common.o
	g++ src/EGamma_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_PIDsep
EEGamma_PIDsep: src/EEGamma_PIDsep.cxx include/EEGamma_PIDsep.h common.o
	g++ src/EEGamma_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_PIDsep
EDopp_PIDsep: src/EDopp_PIDsep.cxx include/EDopp_PIDsep.h common.o
	g++ src/EDopp_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_PIDsep
EEDopp_PIDsep: src/EEDopp_PIDsep.cxx include/EEDopp_PIDsep.h common.o
	g++ src/EEDopp_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_PIDsep
EDopp_ETIPtot_PIDsep: src/EDopp_ETIPtot_PIDsep.cxx include/EDopp_ETIPtot_PIDsep.h common.o
	g++ src/EDopp_ETIPtot_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_ETIPtot_PIDsep
EGamma_PIDsep_mca: src/EGamma_PIDsep_mca.cxx include/EGamma_PIDsep_mca.h common.o
	g++ src/EGamma_PIDsep_mca.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_PIDsep_mca
EDopp_PIDsep_mca: src/EDopp_PIDsep_mca.cxx include/EDopp_PIDsep_mca.h common.o
	g++ src/EDopp_PIDsep_mca.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_PIDsep_mca
TigEE_PIDsep_mca: src/TigEE_PIDsep_mca.cxx include/TigEE_PIDsep_mca.h common.o
	g++ src/TigEE_PIDsep_mca.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TigEE_PIDsep_mca
EGamma_mca_SMOL: src/EGamma_mca_SMOL.cxx include/EGamma_mca_SMOL.h common.o
	g++ src/EGamma_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_mca_SMOL
EGamma_noAB_mca_SMOL: src/EGamma_noAB_mca_SMOL.cxx include/EGamma_noAB_mca_SMOL.h common.o
	g++ src/EGamma_noAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_noAB_mca_SMOL
EGamma_no90deg_mca_SMOL: src/EGamma_no90deg_mca_SMOL.cxx include/EGamma_no90deg_mca_SMOL.h common.o
	g++ src/EGamma_no90deg_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_no90deg_mca_SMOL
EGamma_mca_SMOL_nopos15: src/EGamma_mca_SMOL_nopos15.cxx include/EGamma_mca_SMOL_nopos15.h common.o
	g++ src/EGamma_mca_SMOL_nopos15.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_mca_SMOL_nopos15
EGamma_phi_mca_SMOL: src/EGamma_phi_mca_SMOL.cxx include/EGamma_phi_mca_SMOL.h common.o
	g++ src/EGamma_phi_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_phi_mca_SMOL
EEGamma_mca_SMOL: src/EEGamma_mca_SMOL.cxx include/EEGamma_mca_SMOL.h common.o
	g++ src/EEGamma_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_mca_SMOL
EEGamma_mult2timing_mca_SMOL: src/EEGamma_mult2timing_mca_SMOL.cxx include/EEGamma_mult2timing_mca_SMOL.h common.o
	g++ src/EEGamma_mult2timing_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_mult2timing_mca_SMOL
EEGamma_noAB_mult2timing_mca_SMOL: src/EEGamma_noAB_mult2timing_mca_SMOL.cxx include/EEGamma_noAB_mult2timing_mca_SMOL.h common.o
	g++ src/EEGamma_noAB_mult2timing_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_noAB_mult2timing_mca_SMOL
EEGamma_noAB_noSupp_mult2timing_mca_SMOL: src/EEGamma_noAB_noSupp_mult2timing_mca_SMOL.cxx include/EEGamma_noAB_noSupp_mult2timing_mca_SMOL.h common.o
	g++ src/EEGamma_noAB_noSupp_mult2timing_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_noAB_noSupp_mult2timing_mca_SMOL
EEGamma_ANDgate_mca_SMOL: src/EEGamma_ANDgate_mca_SMOL.cxx include/EEGamma_ANDgate_mca_SMOL.h common.o
	g++ src/EEGamma_ANDgate_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_ANDgate_mca_SMOL
EEGamma_no90deg_mca_SMOL: src/EEGamma_no90deg_mca_SMOL.cxx include/EEGamma_no90deg_mca_SMOL.h common.o
	g++ src/EEGamma_no90deg_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_no90deg_mca_SMOL
EEGamma_only180degcoinc_mca_SMOL: src/EEGamma_only180degcoinc_mca_SMOL.cxx include/EEGamma_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEGamma_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_only180degcoinc_mca_SMOL
EEGamma_only180degcoincnoAB_mca_SMOL: src/EEGamma_only180degcoincnoAB_mca_SMOL.cxx include/EEGamma_only180degcoincnoAB_mca_SMOL.h common.o
	g++ src/EEGamma_only180degcoincnoAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_only180degcoincnoAB_mca_SMOL
EDopp_mca_SMOL: src/EDopp_mca_SMOL.cxx include/EDopp_mca_SMOL.h common.o
	g++ src/EDopp_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_mca_SMOL
EDopp_noAB_mca_SMOL: src/EDopp_noAB_mca_SMOL.cxx include/EDopp_noAB_mca_SMOL.h common.o
	g++ src/EDopp_noAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_noAB_mca_SMOL
EDopp_noAB_mult2up_mca_SMOL: src/EDopp_noAB_mult2up_mca_SMOL.cxx include/EDopp_noAB_mult2up_mca_SMOL.h common.o
	g++ src/EDopp_noAB_mult2up_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_noAB_mult2up_mca_SMOL
EDopp_phi_mca_SMOL: src/EDopp_phi_mca_SMOL.cxx include/EDopp_phi_mca_SMOL.h common.o
	g++ src/EDopp_phi_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_phi_mca_SMOL
EEDopp_mca_SMOL: src/EEDopp_mca_SMOL.cxx include/EEDopp_mca_SMOL.h common.o
	g++ src/EEDopp_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_mca_SMOL
EEDopp_noAB_mca_SMOL: src/EEDopp_noAB_mca_SMOL.cxx include/EEDopp_noAB_mca_SMOL.h common.o
	g++ src/EEDopp_noAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_noAB_mca_SMOL
EGamma_only180degcoinc_mca_SMOL: src/EGamma_only180degcoinc_mca_SMOL.cxx include/EGamma_only180degcoinc_mca_SMOL.h common.o
	g++ src/EGamma_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_only180degcoinc_mca_SMOL
EDopp_only180degcoinc_mca_SMOL: src/EDopp_only180degcoinc_mca_SMOL.cxx include/EDopp_only180degcoinc_mca_SMOL.h common.o
	g++ src/EDopp_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDopp_only180degcoinc_mca_SMOL
EEDopp_only180degcoinc_mca_SMOL: src/EEDopp_only180degcoinc_mca_SMOL.cxx include/EEDopp_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEDopp_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_only180degcoinc_mca_SMOL
EEEGamma_only180degcoinc_mca_SMOL: src/EEEGamma_only180degcoinc_mca_SMOL.cxx include/EEEGamma_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEEGamma_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEEGamma_only180degcoinc_mca_SMOL
EEEDopp_only180degcoinc_mca_SMOL: src/EEEDopp_only180degcoinc_mca_SMOL.cxx include/EEEDopp_only180degcoinc_mca_SMOL.h common.o
	g++ src/EEEDopp_only180degcoinc_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEEDopp_only180degcoinc_mca_SMOL
EEDopp_ANDGate_mca_SMOL: src/EEDopp_ANDGate_mca_SMOL.cxx include/EEDopp_ANDGate_mca_SMOL.h common.o
	g++ src/EEDopp_ANDGate_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_ANDGate_mca_SMOL
EEDopp_ORGate_mca_SMOL: src/EEDopp_ORGate_mca_SMOL.cxx include/EEDopp_ORGate_mca_SMOL.h common.o
	g++ src/EEDopp_ORGate_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_ORGate_mca_SMOL
EEDopp_ORGate_noAB_mca_SMOL: src/EEDopp_ORGate_noAB_mca_SMOL.cxx include/EEDopp_ORGate_noAB_mca_SMOL.h common.o
	g++ src/EEDopp_ORGate_noAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_ORGate_noAB_mca_SMOL
EEDopp_no90deg_mca_SMOL: src/EEDopp_no90deg_mca_SMOL.cxx include/EEDopp_no90deg_mca_SMOL.h common.o
	g++ src/EEDopp_no90deg_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEDopp_no90deg_mca_SMOL
EDoppEGamma_mca_SMOL: src/EDoppEGamma_mca_SMOL.cxx include/EDoppEGamma_mca_SMOL.h common.o
	g++ src/EDoppEGamma_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDoppEGamma_mca_SMOL
EDoppGateAngleEGamma_mca_SMOL: src/EDoppGateAngleEGamma_mca_SMOL.cxx include/EDoppGateAngleEGamma_mca_SMOL.h common.o
	g++ src/EDoppGateAngleEGamma_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EDoppGateAngleEGamma_mca_SMOL
EGammaEDopp_mca_SMOL: src/EGammaEDopp_mca_SMOL.cxx include/EGammaEDopp_mca_SMOL.h common.o
	g++ src/EGammaEDopp_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGammaEDopp_mca_SMOL
EGammaEDopp_noAB_mca_SMOL: src/EGammaEDopp_noAB_mca_SMOL.cxx include/EGammaEDopp_noAB_mca_SMOL.h common.o
	g++ src/EGammaEDopp_noAB_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGammaEDopp_noAB_mca_SMOL
ETIP_mca_SMOL: src/ETIP_mca_SMOL.cxx include/ETIP_mca_SMOL.h common.o
	g++ src/ETIP_mca_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ETIP_mca_SMOL
TTGamma_SMOL: src/TTGamma_SMOL.cxx include/TTGamma_SMOL.h common.o
	g++ src/TTGamma_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TTGamma_SMOL
TDoppTGamma_SMOL: src/TDoppTGamma_SMOL.cxx include/TDoppTGamma_SMOL.h common.o
	g++ src/TDoppTGamma_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TDoppTGamma_SMOL
TTGamma_PIDsep: src/TTGamma_PIDsep.cxx include/TTGamma_PIDsep.h common.o
	g++ src/TTGamma_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TTGamma_PIDsep
TTipTGamma_PIDsep: src/TTipTGamma_PIDsep.cxx include/TTipTGamma_PIDsep.h common.o
	g++ src/TTipTGamma_PIDsep.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TTipTGamma_PIDsep
TIPRingPairs: src/TIPRingPairs.cxx include/TIPRingPairs.h common.o
	g++ src/TIPRingPairs.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o TIPRingPairs
SeparatorTimingPID: src/SeparatorTimingPID.cxx include/SeparatorTimingPID.h include/evt_fmt.h common.o
	g++ src/SeparatorTimingPID.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SeparatorTimingPID
SeparatorSource: src/SeparatorSource.cxx include/SeparatorSource.h include/evt_fmt.h common.o
	g++ src/SeparatorSource.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SeparatorSource
clean:
	rm -rf *~ SortDiagnosticsGRSI SortDiagnosticsSMOL DrawTIPWaveforms PrintRFData ComptonAsym_mca_SMOL ComptonAsym_mca_SMOL_EEGamma ComptonAsym_mca_SMOL_EEDopp ComptonAsym_mca_SMOL_EDopp ComptonAngle_SMOL ComptonAngle_SMOL_EDopp ComptonAngle_SMOL_EGammaEDopp RDCO_SMOL RDCO_SMOL_EDopp RDCO_SMOL_EGammaEDopp AngCorrSMOL AngDistSMOL SumSMOL CheckTimingWindows PlotTimeSepPID EGamma_PIDsep EEGamma_PIDsep EDopp_PIDsep EEDopp_PIDsep EDopp_ETIPtot_PIDsep EGamma_PIDsep_mca EDopp_PIDsep_mca TigEE_PIDsep_mca EGamma_mca_SMOL EGamma_noAB_mca_SMOL EGamma_no90deg_mca_SMOL EGamma_phi_mca_SMOL EEGamma_mca_SMOL EEGamma_mult2timing_mca_SMOL EEGamma_noAB_mult2timing_mca_SMOL EEGamma_noAB_noSupp_mult2timing_mca_SMOL EEGamma_ANDgate_mca_SMOL EEGamma_no90deg_mca_SMOL EEGamma_only180degcoinc_mca_SMOL EEGamma_only180degcoincnoAB_mca_SMOL EDopp_mca_SMOL EDopp_noAB_mca_SMOL EDopp_noAB_mult2up_mca_SMOL EDopp_phi_mca_SMOL EGamma_mca_SMOL_nopos15 EEDopp_mca_SMOL EEDopp_noAB_mca_SMOL EGamma_only180degcoinc_mca_SMOL EDopp_only180degcoinc_mca_SMOL EEDopp_only180degcoinc_mca_SMOL EEEGamma_only180degcoinc_mca_SMOL EEEDopp_only180degcoinc_mca_SMOL EEDopp_ANDGate_mca_SMOL EEDopp_ORGate_mca_SMOL EEDopp_ORGate_noAB_mca_SMOL EEDopp_no90deg_mca_SMOL EDoppEGamma_mca_SMOL EDoppGateAngleEGamma_mca_SMOL EGammaEDopp_mca_SMOL EGammaEDopp_noAB_mca_SMOL ETIP_mca_SMOL TTGamma_SMOL TDoppTGamma_SMOL TTGamma_PIDsep TTipTGamma_PIDsep TIPRingPairs SeparatorTimingPID SeparatorSource common.o
