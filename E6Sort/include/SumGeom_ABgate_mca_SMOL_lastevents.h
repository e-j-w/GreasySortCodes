#ifndef SumGeom_ABgate_mca_SMOL_lastevents_h
#define SumGeom_ABgate_mca_SMOL_lastevents_h

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

double gateELow[MAX_NUM_GATES], gateEHigh[MAX_NUM_GATES];
uint8_t numPctToSort;
double pctToSort[MAX_NUM_PCTTOSORT];
uint64_t sortingSubset; //bit-pattern describing which subsets (%s) of data are currently being sorted 
uint64_t evtsToSort[MAX_NUM_PCTTOSORT], maxEvtsToSort;
uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint32_t numFilesWritten;
uint64_t totalEntriesRead, totalEntriesInFileList;
double addbackE[NGRIFPOS],maxABHitE[NGRIFPOS];
double addbackT[NGRIFPOS];
uint8_t addbackTS[NGRIFPOS];
uint8_t addbackNumCFDFail[NGRIFPOS];

uint8_t prevEvtABPos, currentEvtABPos;

uint64_t numEvtCoinc[MAX_NUM_PCTTOSORT];
uint64_t numEvtCoincSum[MAX_NUM_PCTTOSORT];

double coincGateMin, coincGateMax, coincGate1CFDFailMin, coincGate1CFDFailMax, coincGate2CFDFailMin, coincGate2CFDFailMax;
double sumGateMin, sumGateMax, sumGateCFDMin, sumGateCFDMax, sumGate1CFDFailMin, sumGate1CFDFailMax, sumGate2CFDFailMin, sumGate2CFDFailMax; 
double tRandGateMin, tRandGateMax, leCoincGateMin, leCoincGateMax, leTRandGateMin, leTRandGateMax;

#endif

