CPPFLAGS = -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -Wall

all: GammaTrackingTIGRESS.o GammaTrackingMakeMap GammaTrackingMakeBasis GammaTrackingSortExample GammaTrackingSortFromBasisGridSearch GammaTrackingSortHybrid

GammaTrackingTIGRESS.o: GammaTrackingTIGRESS.cxx GammaTrackingTIGRESS.h
	g++ $(CPPFLAGS) -c -o GammaTrackingTIGRESS.o GammaTrackingTIGRESS.cxx

GammaTrackingMakeMap: gamma_tracking_make_map.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_make_map.cxx $(CPPFLAGS) -o GammaTrackingMakeMap GammaTrackingTIGRESS.o

GammaTrackingMakeBasis: gamma_tracking_make_basis.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_make_basis.cxx $(CPPFLAGS) -o GammaTrackingMakeBasis GammaTrackingTIGRESS.o
	
GammaTrackingSortExample: gamma_tracking_sort_example.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_example.cxx $(CPPFLAGS) -o GammaTrackingSortExample GammaTrackingTIGRESS.o

GammaTrackingSortFromBasisGridSearch: gamma_tracking_sort_from_basis_gridsearch.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_from_basis_gridsearch.cxx $(CPPFLAGS) -o GammaTrackingSortFromBasisGridSearch GammaTrackingTIGRESS.o

GammaTrackingSortHybrid: gamma_tracking_sort_hybrid.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_hybrid.cxx $(CPPFLAGS) -o GammaTrackingSortHybrid GammaTrackingTIGRESS.o
	
clean:
	rm GammaTrackingTIGRESS.o GammaTrackingMakeMap GammaTrackingMakeBasis GammaTrackingSortExample GammaTrackingSortFromBasisGridSearch GammaTrackingSortHybrid
