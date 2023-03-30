#ifndef TDoppTGamma_SMOL_h
#define TDoppTGamma_SMOL_h

#include <iostream>
#include <iomanip>

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

PIDGates *gates;

class TDoppTGamma_SMOL {

	public :

		TDoppTGamma_SMOL(){;} 
		void SortData(const char*,const char*);

};
#endif
