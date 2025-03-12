#ifndef EEGamma_noAB_mca_SMOL_h
#define EEGamma_noAB_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[S32K]; //output .fmca data, sp 0 is sum

class EEGamma_noAB_mca_SMOL{
	public :

		EEGamma_noAB_mca_SMOL(){;}
		void WriteData(const char*);
		uint64_t SortData(const char*, const double, const double, const double);

};

#endif

