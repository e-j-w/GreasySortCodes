#ifndef EEGamma_noAB_mca_SMOL_splitruns_h
#define EEGamma_noAB_mca_SMOL_splitruns_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[9][S32K]; //output .dmca data
// sp 0 is gated energy spectrum
// sp 1 is 180 degree projection
// sp 2 is 180 degree sum
// sp 3 is time random
// sp 4 is time random 180 degree projection
// sp 5 is time random 180 degree sum
// sp 6 is singles (ungated)
// sp 7 is ungated 180 degree projection
// sp 8 is ungated 180 degree sum


class EEGamma_noAB_mca_SMOL_splitruns{
	public :

		EEGamma_noAB_mca_SMOL_splitruns(){;}
		uint64_t SortData(const char*, const char*, const double, const double, const double, const uint8_t, const uint64_t, const double, const double, const double);

};

#endif

