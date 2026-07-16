#ifndef RedistE_Exp_h
#define RedistE_Exp_h

#include <iostream>
#include <iomanip>

using namespace std;

#define NSPECT 100 //maximum number of spectra in an input file
#define MAX_SP_REDIST 8 //maximum number of spectra that can be redistributed

//spectra
double mcaIn[NSPECT][S32K], mcaOut[NSPECT][S32K]; //output .dmca data, sp 0 is sum, sp 1 is 180 degree projection, sp 2 is 180 degree sum 

#endif

