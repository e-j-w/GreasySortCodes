#ifndef EGamma_mca_SMOL_h
#define EGamma_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[S32K]; //output .dmca data, sp 0 is sum

class EGamma_mca_SMOL{
	public :
		EGamma_mca_SMOL(){;}
		void WriteData(const char*);
		uint64_t SortData(const char*, const double);

};

#endif

