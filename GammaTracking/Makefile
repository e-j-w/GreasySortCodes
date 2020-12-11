CPPFLAGS = -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT

all: gamma_tracking_lib gamma_tracking_make_map gamma_tracking_make_basis gamma_tracking_sort_example gamma_tracking_sort_from_basis_gridsearch gamma_tracking_sort_hybrid

gamma_tracking_lib: GammaTrackingTIGRESS.cxx GammaTrackingTIGRESS.h
	g++ $(CPPFLAGS) -c -o libGammaTrackingTIGRESS.o GammaTrackingTIGRESS.cxx

gamma_tracking_make_map: gamma_tracking_make_map.cxx libGammaTrackingTIGRESS.o
	g++ gamma_tracking_make_map.cxx $(CPPFLAGS) -o GammaTrackingMakeMap libGammaTrackingTIGRESS.o

gamma_tracking_make_basis: gamma_tracking_make_basis.cxx libGammaTrackingTIGRESS.o
	g++ gamma_tracking_make_basis.cxx $(CPPFLAGS) -o GammaTrackingMakeBasis libGammaTrackingTIGRESS.o
	
gamma_tracking_sort_example: gamma_tracking_sort_example.cxx
	g++ gamma_tracking_sort_example.cxx $(CPPFLAGS) -o GammaTrackingSortExample libGammaTrackingTIGRESS.o

gamma_tracking_sort_from_basis_gridsearch: gamma_tracking_sort_from_basis_gridsearch.cxx libGammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_from_basis_gridsearch.cxx $(CPPFLAGS) -o GammaTrackingSortFromBasisGridSearch libGammaTrackingTIGRESS.o

gamma_tracking_sort_hybrid: gamma_tracking_sort_hybrid.cxx libGammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_hybrid.cxx $(CPPFLAGS) -o GammaTrackingSortHybrid libGammaTrackingTIGRESS.o
	
clean:
	rm libGammaTrackingTIGRESS.o GammaTrackingMakeMap GammaTrackingMakeBasis GammaTrackingSortExample GammaTrackingSortFromBasisGridSearch GammaTrackingSortHybrid *~
