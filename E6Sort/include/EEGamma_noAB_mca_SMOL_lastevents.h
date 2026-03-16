#ifndef EEGamma_noAB_mca_SMOL_lastevents_h
#define EEGamma_noAB_mca_SMOL_lastevents_h

#include <iostream>
#include <iomanip>

using namespace std;

//spectra
double mcaOut[15][S32K]; //output .dmca data
// sp 0 is gated energy spectrum
// sp 1 is 180 degree projection
// sp 2 is 180 degree sum
// sp 3 is time random
// sp 4 is time random 180 degree projection
// sp 5 is time random 180 degree sum
// sp 6 is singles (ungated)
// sp 7 is ungated 180 degree projection
// sp 8 is ungated 180 degree sum
// sp 9 is leading edge gated energy spectrum
// sp 10 is leading edge 180 degree projection
// sp 11 is leading edge 180 degree sum
// sp 12 is leading edge time random
// sp 13 is leading edge time random 180 degree projection
// sp 14 is leading edge time random 180 degree sum

enum sp_enum{
SP_GATED, SP_SUMOUT, SP_SUMIN,
SP_TR_GATED, SP_TR_SUMOUT, SP_TR_SUMIN,
SP_SINGLES, SP_SINGLES_SUMOUT, SP_SINGLES_SUMIN,
SP_LE_GATED, SP_LE_SUMOUT, SP_LE_SUMIN,
SP_LE_TR_GATED, SP_LE_TR_SUMOUT, SP_LE_TR_SUMIN,
SP_ENUM_LENGTH
};

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint32_t numFilesWritten;
uint64_t totalEntriesRead, totalEntriesInFileList;

#endif

