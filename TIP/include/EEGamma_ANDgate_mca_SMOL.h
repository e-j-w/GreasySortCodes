#ifndef EEGamma_ANDgate_mca_SMOL_h
#define EEGamma_ANDgate_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_GATES 4 //maximum number of energy gates that can be defined

//spectra
float mcaOut[NTIGRING+NTIGSEGRING+1][S32K]; //output .fmca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

class EEGamma_ANDgate_mca_SMOL{
	public :

		EEGamma_ANDgate_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, const double, const double[MAX_GATES], const double[MAX_GATES], const double);

};

#endif

