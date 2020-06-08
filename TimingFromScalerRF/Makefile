all: TimingFromScalerRF

TimingFromScalerRF: TimingFromScalerRF.cxx
	g++ TimingFromScalerRF.cxx -std=c++0x -I${GRSISYS}/include -L${GRSISYS}/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I${GRSISYS}/GRSIData/include -L${GRSISYS}/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o TimingFromScalerRF

clean:
	@echo Cleaning up...
	rm -rf *~ *.o TimingFromScalerRF *tmpdatafile*
