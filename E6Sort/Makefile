LINKFLAGS = -Wl,--no-as-needed $(shell root-config --cflags --libs --glibs) -L/opt/local/lib
CXXFLAGS = -O2 -Wl,--copy-dt-needed-entries
GRSISORT = -L/opt/local/lib -lX11 -lXpm $(shell grsi-config --cflags --all-libs --GRSIData-libs) -I$(GRSISYS)/GRSIData/include
#GRSISORT = -I/data1/jwilliams/Code/GRSISort/include  -L/data1/jwilliams/Code/GRSISort/lib -lTFormat -lTLoops -lTKinematics -lTSRIM -lTBetaDecay -lTCal -lTReaction -lTPeakFitting -lTSuppressed -lTBgo -lTPulseAnalyzer -lTGRSIFit -lTNucleus -lTRawFile -lTDataParser -lTGUI -lXMLParser -lXMLIO -lX11 -lXpm -lGuiHtml -lMinuit -lSpectrum -lMathMore -L/data1/jwilliams/Code/GRSISort/GRSIData/lib -lTAngularCorrelation -lTGRSIDataParser -lTGRSIFormat -lTGenericDetector -lTAries -lTCSM -lTDescant -lTGriffin -lTLaBr -lTMidas -lTPaces -lTRF -lTS3 -lTSceptar -lTSharc -lTSiLi -lTTAC -lTTigress -lTTip -lTTriFoil -lTZeroDegree -lTEmma -lTTrific -lTRcmp -lTSharc2 -I/data1/jwilliams/Code/GRSISort/GRSIData/include -pthread -std=c++17 -m64 -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-variable -I/data1/jwilliams/Code/root_6.28.12/include -L/data1/jwilliams/Code/root_6.28.12/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/data1/jwilliams/Code/root_6.28.12/lib -pthread -lm -ldl -rdynamic -L/usr/lib/x86_64-linux-gnu

all: common.o lin_eq_solver.o AngCorrSMOL EGamma_noAB_mca_SMOL EEGamma_noAB_mca_SMOL EEGamma_noAB_mca_SMOL_splitruns EEGamma_ABgate_mca_SMOL_splitruns EEGamma_noAB_mca_SMOL_lastevents EEGamma_ABgate_mca_SMOL_lastevents SumGeom_ABgate_mca_SMOL_splitruns SumGeom_ABgate_mca_SMOL_lastevents SumGeom_noAB_mca_SMOL_lastevents EEGamma_AB_mca_SMOL_lastevents EEGamma_noAB_sumEgate_mca_SMOL EEGamma_modAB_sumEgate_mca_SMOL EEGamma_modAB_comptonCorr_mca_SMOL EEGamma_modAB_mca_SMOL EGamma_modAB_mca_SMOL SortDiagnosticsSMOL DecayCurveSMOL GainCurveSMOL SeparatorSource ReprocessGain_SMOL RemoveCore_SMOL SplitData_SMOL RedistE_Exp

common.o: src/common.cxx src/position_arrays.cxx include/common.h include/evt_fmt.h
	g++ src/common.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -c -o common.o
lin_eq_solver.o: src/lin_eq_solver.c include/lin_eq_solver.h
	g++ src/lin_eq_solver.c -I./include $(CXXFLAGS) -c -o lin_eq_solver.o
AngCorrSMOL: src/AngCorrSMOL.cxx include/AngCorrSMOL.h include/evt_fmt.h common.o
	g++ src/AngCorrSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o AngCorrSMOL
EGamma_noAB_mca_SMOL: src/EGamma_noAB_mca_SMOL.cxx include/EGamma_noAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EGamma_noAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EGamma_noAB_mca_SMOL
EEGamma_noAB_mca_SMOL: src/EEGamma_noAB_mca_SMOL.cxx include/EEGamma_noAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_noAB_mca_SMOL
EEGamma_noAB_mca_SMOL_splitruns: src/EEGamma_noAB_mca_SMOL_splitruns.cxx include/EEGamma_noAB_mca_SMOL_splitruns.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_mca_SMOL_splitruns.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_noAB_mca_SMOL_splitruns
EEGamma_ABgate_mca_SMOL_splitruns: src/EEGamma_ABgate_mca_SMOL_splitruns.cxx include/EEGamma_ABgate_mca_SMOL_splitruns.h include/evt_fmt.h common.o
	g++ src/EEGamma_ABgate_mca_SMOL_splitruns.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_ABgate_mca_SMOL_splitruns
EEGamma_noAB_mca_SMOL_lastevents: src/EEGamma_noAB_mca_SMOL_lastevents.cxx include/EEGamma_noAB_mca_SMOL_lastevents.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_mca_SMOL_lastevents.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_noAB_mca_SMOL_lastevents
EEGamma_ABgate_mca_SMOL_lastevents: src/EEGamma_ABgate_mca_SMOL_lastevents.cxx include/EEGamma_ABgate_mca_SMOL_lastevents.h include/evt_fmt.h common.o
	g++ src/EEGamma_ABgate_mca_SMOL_lastevents.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_ABgate_mca_SMOL_lastevents
SumGeom_ABgate_mca_SMOL_splitruns: src/SumGeom_ABgate_mca_SMOL_splitruns.cxx include/SumGeom_ABgate_mca_SMOL_splitruns.h include/evt_fmt.h common.o
	g++ src/SumGeom_ABgate_mca_SMOL_splitruns.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o SumGeom_ABgate_mca_SMOL_splitruns
SumGeom_ABgate_mca_SMOL_lastevents: src/SumGeom_ABgate_mca_SMOL_lastevents.cxx include/SumGeom_ABgate_mca_SMOL_lastevents.h include/evt_fmt.h common.o
	g++ src/SumGeom_ABgate_mca_SMOL_lastevents.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o SumGeom_ABgate_mca_SMOL_lastevents
SumGeom_noAB_mca_SMOL_lastevents: src/SumGeom_noAB_mca_SMOL_lastevents.cxx include/SumGeom_noAB_mca_SMOL_lastevents.h include/evt_fmt.h common.o
	g++ src/SumGeom_noAB_mca_SMOL_lastevents.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o SumGeom_noAB_mca_SMOL_lastevents
EEGamma_AB_mca_SMOL_lastevents: src/EEGamma_AB_mca_SMOL_lastevents.cxx include/EEGamma_AB_mca_SMOL_lastevents.h include/evt_fmt.h common.o
	g++ src/EEGamma_AB_mca_SMOL_lastevents.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_AB_mca_SMOL_lastevents
EEGamma_noAB_sumEgate_mca_SMOL: src/EEGamma_noAB_sumEgate_mca_SMOL.cxx include/EEGamma_noAB_sumEgate_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_sumEgate_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_noAB_sumEgate_mca_SMOL
EEGamma_modAB_sumEgate_mca_SMOL: src/EEGamma_modAB_sumEgate_mca_SMOL.cxx include/EEGamma_modAB_sumEgate_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_modAB_sumEgate_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_modAB_sumEgate_mca_SMOL
EEGamma_modAB_comptonCorr_mca_SMOL: src/EEGamma_modAB_comptonCorr_mca_SMOL.cxx include/EEGamma_modAB_comptonCorr_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_modAB_comptonCorr_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_modAB_comptonCorr_mca_SMOL
EEGamma_modAB_mca_SMOL: src/EEGamma_modAB_mca_SMOL.cxx include/EEGamma_modAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_modAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EEGamma_modAB_mca_SMOL
EGamma_modAB_mca_SMOL: src/EGamma_modAB_mca_SMOL.cxx include/EGamma_modAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EGamma_modAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o EGamma_modAB_mca_SMOL
SortDiagnosticsSMOL: src/SortDiagnosticsSMOL.cxx include/SortDiagnosticsSMOL.h include/evt_fmt.h common.o
	g++ src/SortDiagnosticsSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o SortDiagnosticsSMOL
DecayCurveSMOL: src/DecayCurveSMOL.cxx include/DecayCurveSMOL.h include/evt_fmt.h common.o
	g++ src/DecayCurveSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o DecayCurveSMOL
GainCurveSMOL: src/GainCurveSMOL.cxx include/GainCurveSMOL.h include/evt_fmt.h common.o
	g++ src/GainCurveSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o GainCurveSMOL
SeparatorSource: src/SeparatorSource.cxx include/SeparatorSource.h include/evt_fmt.h common.o
	g++ src/SeparatorSource.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SeparatorSource
ReprocessGain_SMOL: src/ReprocessGain_SMOL.cxx include/ReprocessGain_SMOL.h include/evt_fmt.h common.o lin_eq_solver.o
	g++ src/ReprocessGain_SMOL.cxx -I./include common.o lin_eq_solver.o $(LINKFLAGS) $(CXXFLAGS) -o ReprocessGain_SMOL
RemoveCore_SMOL: src/RemoveCore_SMOL.cxx include/RemoveCore_SMOL.h include/evt_fmt.h common.o lin_eq_solver.o
	g++ src/RemoveCore_SMOL.cxx -I./include common.o lin_eq_solver.o $(LINKFLAGS) $(CXXFLAGS) -o RemoveCore_SMOL
SplitData_SMOL: src/SplitData_SMOL.cxx include/SplitData_SMOL.h include/evt_fmt.h common.o lin_eq_solver.o
	g++ src/SplitData_SMOL.cxx -I./include common.o lin_eq_solver.o $(LINKFLAGS) $(CXXFLAGS) -o SplitData_SMOL
RedistE_Exp: src/RedistE_Exp.cxx include/RedistE_Exp.h common.o
	g++ src/RedistE_Exp.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) -o RedistE_Exp
clean:
	rm -rf *~ EGamma_noAB_mca_SMOL AngCorrSMOL EEGamma_noAB_mca_SMOL EEGamma_noAB_mca_SMOL_splitruns EEGamma_ABgate_mca_SMOL_splitruns EEGamma_noAB_mca_SMOL_lastevents EEGamma_ABgate_mca_SMOL_lastevents SumGeom_ABgate_mca_SMOL_splitruns SumGeom_ABgate_mca_SMOL_lastevents SumGeom_noAB_mca_SMOL_lastevents EEGamma_AB_mca_SMOL_lastevents EEGamma_noAB_sumEgate_mca_SMOL EEGamma_modAB_sumEgate_mca_SMOL EEGamma_modAB_comptonCorr_mca_SMOL EEGamma_modAB_mca_SMOL EGamma_modAB_mca_SMOL SortDiagnosticsSMOL DecayCurveSMOL GainCurveSMOL SeparatorSource ReprocessGain_SMOL RemoveCore_SMOL SplitData_SMOL RedistE_Exp common.o lin_eq_solver.o
