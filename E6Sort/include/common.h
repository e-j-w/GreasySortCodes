#ifndef common_h
#define common_h

#include <iostream>
#include <iomanip>
#include <string.h>
#include "TVector3.h"

#include "evt_fmt.h"

using namespace std;

#define NHPGERING 6 //number of TIGRESS rings
#define NTIGSEGRING 13 //number of TIGRESS segment rings (last ring contains hits with no segments)
#define NTIGPOS 16 //number of TIGRESS positions (clovers)

#define MAX_NUM_PARTICLE 4 //maximum particle multiplicity

#define MAXNUMHPGEHIT 64 //maximum number of HPGe hits per event

#define AMU 931.4941 //atomic mass unit, MeV/c^2

#define S32K 32768 //maximum number of channels per spectrum in .mca and .fmca (changing breaks file compatibility)

#define MIN_HPGE_EAB 5 //minimum HPGe energy (in keV) for addback

#define PI 3.14159265359

#define ADDBACK_TIMING_GATE 200.0 //maximum time (in ns) that can separate events which will be summed in addback
#define COINC_TIMING_GATE   200.0 //maximum time (in ns) that can separate sorted coincidences (following addback)
#define SUM_TIMING_GATE     200.0 //maximum time (in ns) that can separate sorted coincidences (following addback), affects summing histograms

//GLOBAL VARIABLES
//(static to avoid multiple declaration when linking)

//timing windows
//static Double_t tigtigTGate[2] = {-15, 5}; // super narrow TIGRESS - TIGRESS timing window (ns)
static Double_t hpgehpgeTGate[2] = {-100, 10}; // narrow HPGe - HPGe timing window (ns)
static Double_t hpgehpgeABGate[2] = {-200, 50}; // addback HPGe - HPGe timing window (ns)
static Double_t hpgehpgeTRandGate[2] = {-2000, -1300}; // time-random HPGe - HPGe timing window (ns)
static Double_t tigBGOTGate[2] = {0, 380}; // TIGRESS - BGO timing window (ns)

//FUNCTION PROTOTYPES
TVector3 getGeVector(const uint8_t core, const uint8_t seg, const uint8_t forwardPos);
Double_t getGeHitDistance(const uint8_t core1, const uint8_t seg1, const uint8_t core2, const uint8_t seg2, const uint8_t forwardPos);
bool gate1D(const Double_t value, const Double_t min, const Double_t max);
Int_t getHPGeRing(const float theta);
Int_t getHPGeSegmentRing(const float theta);
Int_t getHPGePhiRing(const float theta, const float phi);

#endif
