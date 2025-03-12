#ifndef EGamma_modAB_mca_SMOL_h
#define EGamma_modAB_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[3][S32K]; //output .dmca data, sp 0 is gated energy spectrum, sp 1 is 180 degree projection, sp 2 is 180 degree sum

class EGamma_modAB_mca_SMOL{
	public :

		EGamma_modAB_mca_SMOL(){;}
		void WriteData(const char*);
		uint64_t SortData(const char*, const double, const double);

};

#endif

