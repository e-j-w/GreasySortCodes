#ifndef common_h
#define common_h

#include <iostream>
#include <iomanip>
#include <string.h>
#include "TVector3.h"

#include "evt_fmt.h"

using namespace std;

#define NTIP 128 //number of TIP channels
#define NTIPRING 10 //number of TIP rings
#define NTIGRING 6 //number of TIGRESS rings
#define NTIGSEGRING 13 //number of TIGRESS segment rings (last ring contains hits with no segments)
#define NTIGPOS 16 //number of TIGRESS positions (clovers)

#define MAX_NUM_PARTICLE 4 //maximum particle multiplicity

#define MAXNUMTIPHIT 30 //maximum number of TIP hits per event
#define MAXNUMTIGHIT 30 //maximum number of TIGRESS hits per event

#define AMU 931.4941 //atomic mass unit, MeV/c^2

#define S32K 32768 //maximum number of channels per spectrum in .mca and .fmca (changing breaks file compatibility)

#define TIPTIPFLAG 61
#define TIGTIGFLAG 62
#define TIPTIGFLAG 63

#define MIN_TIG_EAB 15 //minimum TIGRESS energy (in keV) for addback

#define PI 3.14159265359

//GLOBAL VARIABLES
//(static to avoid multiple declaration when linking)

//timing windows
//static Double_t tigtigTGate[2] = {-15, 5}; // super narrow TIGRESS - TIGRESS timing window (ns)
static Double_t tigtigTGate[2] = {-20, 20}; // narrow TIGRESS - TIGRESS timing window (ns)
static Double_t tigBGOTGate[2] = {0, 380}; // TIGRESS - BGO timing window (ns)

//FUNCTION PROTOTYPES
TVector3 getTigVector(uint8_t core, uint8_t seg);
bool gate1D(const Double_t value, const Double_t min, const Double_t max);
Int_t getTIGRESSRing(const float theta);
Int_t getTIGRESSSegmentRing(const float theta);
Int_t getTIGRESSPhiRing(const float theta, const float phi);

#endif
