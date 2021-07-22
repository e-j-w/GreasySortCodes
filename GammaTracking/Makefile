CPPFLAGS = -std=c++0x -I$(GRSISYS)/include -L$(GRSISYS)/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$(GRSISYS)/GRSIData/include -L$(GRSISYS)/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT

all: gamma_tracking_lib gamma_tracking_make_map gamma_tracking_make_basis gamma_tracking_sort_example gamma_tracking_sort_from_basis_gridsearch gamma_tracking_sort_hybrid

gamma_tracking_lib: GammaTrackingTIGRESS.cxx GammaTrackingTIGRESS.h
	g++ $(CPPFLAGS) -c -o GammaTrackingTIGRESS.o GammaTrackingTIGRESS.cxx

gamma_tracking_make_map: gamma_tracking_make_map.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_make_map.cxx $(CPPFLAGS) -o GammaTrackingMakeMap GammaTrackingTIGRESS.o

gamma_tracking_make_basis: gamma_tracking_make_basis.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_make_basis.cxx $(CPPFLAGS) -o GammaTrackingMakeBasis GammaTrackingTIGRESS.o
	
gamma_tracking_sort_example: gamma_tracking_sort_example.cxx
	g++ gamma_tracking_sort_example.cxx $(CPPFLAGS) -o GammaTrackingSortExample GammaTrackingTIGRESS.o

gamma_tracking_sort_from_basis_gridsearch: gamma_tracking_sort_from_basis_gridsearch.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_from_basis_gridsearch.cxx $(CPPFLAGS) -o GammaTrackingSortFromBasisGridSearch GammaTrackingTIGRESS.o

gamma_tracking_sort_hybrid: gamma_tracking_sort_hybrid.cxx GammaTrackingTIGRESS.o
	g++ gamma_tracking_sort_hybrid.cxx $(CPPFLAGS) -o GammaTrackingSortHybrid GammaTrackingTIGRESS.o
	
clean:
	rm GammaTrackingTIGRESS.o GammaTrackingMakeMap GammaTrackingMakeBasis GammaTrackingSortExample GammaTrackingSortFromBasisGridSearch GammaTrackingSortHybrid *~
