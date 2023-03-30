#ifndef ComptonAngle_S_h
#define ComptonAngle_S_h

#include <iostream>
#include <iomanip>

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

#include "TGraphErrors.h"

using namespace std;

class ComptonAngle_S {

	public :

		ComptonAngle_S(){;} 
		void SortData(const char*);
    void WriteData(const char*);

};
#endif
