all: TimingFromScalerRF

TimingFromScalerRF: TimingFromScalerRF.cxx
	g++ TimingFromScalerRF.cxx -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -o TimingFromScalerRF

clean:
	@echo Cleaning up...
	rm -rf *~ *.o TimingFromScalerRF *tmpdatafile*
