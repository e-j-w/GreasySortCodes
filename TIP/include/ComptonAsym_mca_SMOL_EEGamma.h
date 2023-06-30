#ifndef ComptonAngle_SEEEE_h
#define ComptonAngle_SEEEE_h

#include <iostream>
#include <iomanip>

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

#include "TGraphErrors.h"

using namespace std;

class ComptonAngle_SEE {

	public :

		ComptonAngle_SEE(){;} 
		void SortData(const char*, const double, const double);
    void WriteData(const char*);

};
#endif
