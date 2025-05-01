#ifndef EEGamma_modAB_mca_SMOL_h
#define EEGamma_modAB_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[6][S32K]; //output .dmca data
// sp 0 is gated energy spectrum
// sp 1 is 180 degree projection
// sp 2 is 180 degree sum
// sp 3 is time random
// sp 4 is time random 180 degree projection
// sp 5 is time random 180 degree sum

class EEGamma_modAB_mca_SMOL{
	public :
		EEGamma_modAB_mca_SMOL(){;}
		void WriteData(const char*);
		uint64_t SortData(const char*, const double, const double, const double, const double, const double);

};

#endif

