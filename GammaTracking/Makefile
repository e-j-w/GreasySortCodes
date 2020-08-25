CPPFLAGS = -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT

all: gamma_tracking_make_map gamma_tracking_sort_example
	
gamma_tracking_make_map: gamma_tracking_make_map.cxx
	g++ gamma_tracking_make_map.cxx $(CPPFLAGS) -o GammaTrackingMakeMap
	
gamma_tracking_sort_example: gamma_tracking_sort_example.cxx
	g++ gamma_tracking_sort_example.cxx $(CPPFLAGS) -o GammaTrackingSortExample
	
clean:
	rm GammaTrackingMakeMap GammaTrackingSortExample *~
