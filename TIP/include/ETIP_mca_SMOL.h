#ifndef ETIP_mca_SMOL_h
#define ETIP_mca_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
float mcaOut[NTIPRING+1][S32K]; //output .mca data, sp 0 is sum, sp 1-6 are TIGRESS rings, sp 7-18 are segment rings

class ETIP_mca_SMOL{
	public :

		ETIP_mca_SMOL(){;}
		void WriteData(const char*);
		void SortData(const char*, double);

};

#endif

