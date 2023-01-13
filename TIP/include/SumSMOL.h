#ifndef SumG_h
#define SumG_h

#include <iostream>
#include <iomanip>

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

class Sum {

	public :

		Sum(){;}
		uint64_t SortData(const char *, FILE *);
};
#endif