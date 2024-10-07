#ifndef EEEDopp_only180degcoinc_mca_SMOL_h
#define EEEDopp_only180degcoinc_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

PIDGates *gates;

class EEEDopp_only180degcoinc_mca_SMOL{
	public :

		EEEDopp_only180degcoinc_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, const double, const double, const double);

};

#endif

