#ifndef SumGeom_ABgate_mca_SMOL_splitruns_h
#define SumGeom_ABgate_mca_SMOL_splitruns_h

#include <iostream>
#include <iomanip>

using namespace std;

#define MAX_NUM_GATES 10 //maximum number of energy gates to sort
#define MAX_NUM_SUBSETS 30 //maximum number of data subsets for sorting

#define EVT_MIX_SEARCH_DEPTH 1000 //number of events to search for 'coincidences' in
#define E_THRESHOLD 250 //in keV

uint32_t numSubsets;
uint32_t sortingSubset; //which subset is currently being sorted
uint64_t subsetEvtsToSort;
uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint64_t totalEntriesRead, totalEntriesInFileList;
double addbackE[NGRIFPOS],maxABHitE[NGRIFPOS];
double addbackT[NGRIFPOS];
uint8_t addbackTS[NGRIFPOS];
uint8_t addbackNumCFDFail[NGRIFPOS];

uint64_t gateABHitEvtNum, gateSumHit1EvtNum;
uint8_t prevEvtGateABPos, prevEvtGatedSumHit1Pos, prevEvtGatedSumHit2Pos;
uint64_t sumHit1EvtNum;
uint8_t prevEvtSumHit1Pos, prevEvtSumHit2Pos;

//counters
uint64_t numEvtGatedCoinc[MAX_NUM_SUBSETS];
uint64_t numEvtGatedCoincSum[MAX_NUM_SUBSETS];
uint64_t numEvtCoinc[MAX_NUM_SUBSETS];
uint64_t numEvtCoincSum[MAX_NUM_SUBSETS];
uint64_t numSinglesHits, numCoincPairs, num180DegCoincPairs;


#endif

