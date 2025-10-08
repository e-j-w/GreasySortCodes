#ifndef ReprocessGain_SMOL_h
#define ReprocessGain_SMOL_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_EVT_WINDOWS_PER_TREE  1024
#define MAX_INPUT_E                8
#define EN_WINDOW_WIDTH            5.0 //width of windows used for energy evaluation, in keV
#define MAX_RUN_TIME_SEC           29800.0


class ReprocessGain_SMOL{
	public :

		ReprocessGain_SMOL(){;} 
		void SortData(const char *sfile, const char *efile, const char *outfile, const uint64_t evalWindowSize, const uint64_t tWindowSize, const double en[MAX_INPUT_E], const int numEnVals);
};
#endif
