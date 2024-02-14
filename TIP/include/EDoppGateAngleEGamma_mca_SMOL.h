#ifndef EDoppGateAngleEGamma_mca_SMOL_h
#define EDoppGateAngleEGamma_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

PIDGates *gates;

class EDoppGateAngleEGamma_mca_SMOL{
	public :

		EDoppGateAngleEGamma_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, const double, const double, const double, const double, const double);

};

#endif

