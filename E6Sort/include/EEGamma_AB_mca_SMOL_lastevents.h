#ifndef EEGamma_AB_mca_SMOL_lastevents_h
#define EEGamma_AB_mca_SMOL_lastevents_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_NUM_GATES 10 //maximum number of energy gates to sort
#define MAX_NUM_PCTTOSORT 30 //maximum number of percentage values for sorting

enum sp_enum{
SP_GATED, SP_SUMOUT, SP_SUMIN, SP_SUMOUT_CFD, SP_SUMIN_CFD,
SP_TR_GATED, SP_TR_SUMOUT, SP_TR_SUMIN,
SP_SINGLES, SP_SINGLES_SUMOUT, SP_SINGLES_SUMIN,
SP_LE_GATED, SP_LE_SUMOUT, SP_LE_SUMIN,
SP_LE_TR_GATED, SP_LE_TR_SUMOUT, SP_LE_TR_SUMIN,
SP_ENUM_LENGTH
};

//spectra
double mcaOut[MAX_NUM_GATES][MAX_NUM_PCTTOSORT][SP_ENUM_LENGTH][S32K]; //output .dmca data
// sp 0 is gated energy spectrum
// sp 1 is 180 degree projection
// sp 2 is 180 degree sum
// sp 3 is 180 degree projection using CFD timing for 180 degree hits
// sp 4 is 180 degree sum using CFD timing for 180 degree hits
// sp 5 is time random
// sp 6 is time random 180 degree projection
// sp 7 is time random 180 degree sum
// sp 8 is singles (ungated)
// sp 9 is ungated 180 degree projection
// sp 10 is ungated 180 degree sum
// sp 11 is leading edge gated energy spectrum
// sp 12 is leading edge 180 degree projection
// sp 13 is leading edge 180 degree sum
// sp 14 is leading edge time random
// sp 15 is leading edge time random 180 degree projection
// sp 16 is leading edge time random 180 degree sum


uint8_t numEGates;
double gateELow[MAX_NUM_GATES], gateEHigh[MAX_NUM_GATES];
uint8_t numPctToSort;
double pctToSort[MAX_NUM_PCTTOSORT];
uint64_t sortingSubset; //bit-pattern describing which subsets (%s) of data are currently being sorted 
uint64_t evtsToSort[MAX_NUM_PCTTOSORT], maxEvtsToSort;
uint8_t hitMap180degAB[NGRIFPOS][NGRIFPOS]; //1st index = clover of hit, 2nd index = clover of 2nd hit, val = 1 indicates 180 degree summing occurs
uint32_t numFilesWritten;
uint64_t totalEntriesRead, totalEntriesInFileList;
double addbackE[NGRIFPOS],maxABHitE[NGRIFPOS];
double addbackT[NGRIFPOS];
uint8_t addbackTS[NGRIFPOS];
uint8_t addbackNumCFDFail[NGRIFPOS];

double coincGateMin, coincGateMax, coincGate1CFDFailMin, coincGate1CFDFailMax, coincGate2CFDFailMin, coincGate2CFDFailMax;
double sumGateMin, sumGateMax, sumGateCFDMin, sumGateCFDMax, sumGate1CFDFailMin, sumGate1CFDFailMax, sumGate2CFDFailMin, sumGate2CFDFailMax; 
double tRandGateMin, tRandGateMax, leCoincGateMin, leCoincGateMax, leTRandGateMin, leTRandGateMax;

#endif

