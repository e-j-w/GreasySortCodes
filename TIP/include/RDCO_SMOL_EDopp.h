#ifndef RDCO_S_h
#define RDCO_S_h

#include <iostream>
#include <iomanip>

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

PIDGates *gates;

class RDCO_SMOL_EDopp {

	public :

		RDCO_SMOL_EDopp(){;} 
		void SortData(const char*);

};
#endif
