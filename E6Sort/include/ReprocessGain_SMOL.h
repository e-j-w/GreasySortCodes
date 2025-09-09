#ifndef ReprocessGain_SMOL_h
#define ReprocessGain_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_TIME_WINDOWS_PER_TREE  1024
#define MAX_INPUT_E                8
#define EN_WINDOW_WIDTH            5.0 //width of windows used for energy evaluation, in keV


class ReprocessGain_SMOL{
	public :

		ReprocessGain_SMOL(){;} 
		void SortData(const char *sfile, const char *efile, const char *outfile, const double evalWindowSize, const double tWindowSize, const double en[MAX_INPUT_E], const int numEnVals);
};
#endif
