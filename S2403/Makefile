all: SortCode

SortCode: SortCode.cxx SortCode.h
	g++ SortCode.cxx -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib -I/data4/S2403/GRSISort/include -L/data4/S2403/GRSISort/lib -lTGRSIint -lTFormat -lTGRSIint -lTFormat -lTLoops -lTKinematics -lTSRIM -lTBetaDecay -lTCal -lTReaction -lTPeakFitting -lTSuppressed -lTBgo -lTPulseAnalyzer -lGROOT -lTGRSIFit -lTNucleus -lTRawFile -lTDataParser -lTGUI -lXMLParser -lXMLIO -lX11 -lXpm -lProof -lGuiHtml -lMinuit -lSpectrum -L/data4/S2403/GRSISort/GRSIData/lib -lTAngularCorrelation -lTCSM -lTDescant -lTGRSIDataParser -lTGRSIFormat -lTGenericDetector -lTGriffin -lTLaBr -lTMidas -lTPaces -lTRF -lTS3 -lTSceptar -lTSharc -lTSiLi -lTTAC -lTTigress -lTTip -lTTriFoil -lTZeroDegree -lTEmma -lTTrific -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lSpectrum -lMinuit -o SortData

clean:
	rm -rf *~ SortData
