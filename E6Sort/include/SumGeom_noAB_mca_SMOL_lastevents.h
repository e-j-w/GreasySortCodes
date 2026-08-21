#ifndef SumGeom_noAB_mca_SMOL_lastevents_h
#define SumGeom_noAB_mca_SMOL_lastevents_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_NUM_GATES 10 //maximum number of energy gates to sort
#define MAX_NUM_PCTTOSORT 30 //maximum number of percentage values for sorting

#define EVT_MIX_SEARCH_DEPTH 1000 //number of events to search for 'coincidences' in
#define E_THRESHOLD 250 //in keV

uint8_t numPctToSort;
double pctToSort[MAX_NUM_PCTTOSORT];
uint64_t sortingSubset; //bit-pattern describing which subsets (%s) of data are currently being sorted 
uint64_t evtsToSort[MAX_NUM_PCTTOSORT], maxEvtsToSort;
uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint64_t totalEntriesRead, totalEntriesInFileList;

uint64_t gateEvtNum, sumHit1EvtNum;
uint8_t prevEvtGatePos, prevEvtSumHit1Pos, prevEvtSumHit2Pos;

//counters
uint64_t numEvtCoinc[MAX_NUM_PCTTOSORT];
uint64_t numEvtCoincSum[MAX_NUM_PCTTOSORT];
uint64_t numSinglesHits, numCoincPairs, num180DegCoincPairs;


#endif

