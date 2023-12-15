#ifndef EEDopp_mca_SMOL_ANDGate_h
#define EEDopp_mca_SMOL_ANDGate_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

PIDGates *gates;

class EEDopp_mca_SMOL_ANDGate{
	public :

		EEDopp_mca_SMOL_ANDGate(){;}
		void WriteData(const char*);
		void SortData(char const *sfile, const uint8_t numEGates, const double eLow[10], const double eHigh[10], const double keVPerBin);

};

#endif

