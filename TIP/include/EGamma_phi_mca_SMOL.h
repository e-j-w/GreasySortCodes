#ifndef EGamma_phi_mca_SMOL_h
#define EGamma_phi_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[5][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

class EGamma_phi_mca_SMOL{
	public :

		EGamma_phi_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, double);

};

#endif

