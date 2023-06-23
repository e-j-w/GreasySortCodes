#ifndef EGamma_no90deg_mca_SMOL_h
#define EGamma_no90deg_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

class EGamma_mca_SMOL{
	public :

		EGamma_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, double);

};

#endif

