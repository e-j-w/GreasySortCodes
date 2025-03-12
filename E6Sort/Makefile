LINKFLAGS = -Wl,--no-as-needed $(shell root-config --cflags --libs --glibs) -lSpectrum -lMinuit -lGuiHtml -lTreePlayer -lTMVA -L/opt/local/lib -lX11 -lXpm
CXXFLAGS = -O2 -Wl,--copy-dt-needed-entries
GRSISORT = -L/opt/local/lib -lX11 -lXpm $(shell grsi-config --cflags --all-libs --GRSIData-libs) -I$(GRSISYS)/GRSIData/include
#GRSISORT = -I/data1/jwilliams/Code/GRSISort/include  -L/data1/jwilliams/Code/GRSISort/lib -lTFormat -lTLoops -lTKinematics -lTSRIM -lTBetaDecay -lTCal -lTReaction -lTPeakFitting -lTSuppressed -lTBgo -lTPulseAnalyzer -lTGRSIFit -lTNucleus -lTRawFile -lTDataParser -lTGUI -lXMLParser -lXMLIO -lX11 -lXpm -lGuiHtml -lMinuit -lSpectrum -lMathMore -L/data1/jwilliams/Code/GRSISort/GRSIData/lib -lTAngularCorrelation -lTGRSIDataParser -lTGRSIFormat -lTGenericDetector -lTAries -lTCSM -lTDescant -lTGriffin -lTLaBr -lTMidas -lTPaces -lTRF -lTS3 -lTSceptar -lTSharc -lTSiLi -lTTAC -lTTigress -lTTip -lTTriFoil -lTZeroDegree -lTEmma -lTTrific -lTRcmp -lTSharc2 -I/data1/jwilliams/Code/GRSISort/GRSIData/include -pthread -std=c++17 -m64 -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-variable -I/data1/jwilliams/Code/root_6.28.12/include -L/data1/jwilliams/Code/root_6.28.12/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,/data1/jwilliams/Code/root_6.28.12/lib -pthread -lm -ldl -rdynamic -L/usr/lib/x86_64-linux-gnu

all: common.o lin_eq_solver.o EGamma_mca_SMOL EGamma_noAB_mca_SMOL EEGamma_noAB_mca_SMOL EEGamma_noAB_tRand_mca_SMOL EEGamma_mca_SMOL EEGamma_modAB_mca_SMOL EGamma_modAB_mca_SMOL SortDiagnosticsSMOL DecayCurveSMOL GainCurveSMOL SeparatorSource ReprocessAB_SMOL ReprocessGain_SMOL

common.o: src/common.cxx src/position_arrays.cxx include/common.h include/evt_fmt.h
	g++ src/common.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -c -o common.o
lin_eq_solver.o: src/lin_eq_solver.c include/lin_eq_solver.h
	g++ src/lin_eq_solver.c -I./include $(CXXFLAGS) -c -o lin_eq_solver.o
EGamma_mca_SMOL: src/EGamma_mca_SMOL.cxx include/EGamma_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EGamma_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_mca_SMOL
EGamma_noAB_mca_SMOL: src/EGamma_noAB_mca_SMOL.cxx include/EGamma_noAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EGamma_noAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_noAB_mca_SMOL
EEGamma_noAB_mca_SMOL: src/EEGamma_noAB_mca_SMOL.cxx include/EEGamma_noAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_noAB_mca_SMOL
EEGamma_noAB_tRand_mca_SMOL: src/EEGamma_noAB_tRand_mca_SMOL.cxx include/EEGamma_noAB_tRand_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_noAB_tRand_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_noAB_tRand_mca_SMOL
EEGamma_mca_SMOL: src/EEGamma_mca_SMOL.cxx include/EEGamma_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_mca_SMOL
EEGamma_modAB_mca_SMOL: src/EEGamma_modAB_mca_SMOL.cxx include/EEGamma_modAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EEGamma_modAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EEGamma_modAB_mca_SMOL
EGamma_modAB_mca_SMOL: src/EGamma_modAB_mca_SMOL.cxx include/EGamma_modAB_mca_SMOL.h include/evt_fmt.h common.o
	g++ src/EGamma_modAB_mca_SMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o EGamma_modAB_mca_SMOL
SortDiagnosticsSMOL: src/SortDiagnosticsSMOL.cxx include/SortDiagnosticsSMOL.h include/evt_fmt.h common.o
	g++ src/SortDiagnosticsSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SortDiagnosticsSMOL
DecayCurveSMOL: src/DecayCurveSMOL.cxx include/DecayCurveSMOL.h include/evt_fmt.h common.o
	g++ src/DecayCurveSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o DecayCurveSMOL
GainCurveSMOL: src/GainCurveSMOL.cxx include/GainCurveSMOL.h include/evt_fmt.h common.o
	g++ src/GainCurveSMOL.cxx -I./include $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o GainCurveSMOL
SeparatorSource: src/SeparatorSource.cxx include/SeparatorSource.h include/evt_fmt.h common.o
	g++ src/SeparatorSource.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o SeparatorSource
ReprocessAB_SMOL: src/ReprocessAB_SMOL.cxx include/ReprocessAB_SMOL.h include/evt_fmt.h common.o
	g++ src/ReprocessAB_SMOL.cxx -I./include common.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ReprocessAB_SMOL
ReprocessGain_SMOL: src/ReprocessGain_SMOL.cxx include/ReprocessGain_SMOL.h include/evt_fmt.h common.o lin_eq_solver.o
	g++ src/ReprocessGain_SMOL.cxx -I./include common.o lin_eq_solver.o $(LINKFLAGS) $(CXXFLAGS) $(GRSISORT) -o ReprocessGain_SMOL
clean:
	rm -rf *~ EGamma_mca_SMOL EGamma_noAB_mca_SMOL EEGamma_noAB_mca_SMOL EEGamma_noAB_tRand_mca_SMOL EEGamma_mca_SMOL EEGamma_modAB_mca_SMOL EGamma_modAB_mca_SMOL SortDiagnosticsSMOL DecayCurveSMOL GainCurveSMOL SeparatorSource ReprocessAB_SMOL ReprocessGain_SMOL common.o lin_eq_solver.o
