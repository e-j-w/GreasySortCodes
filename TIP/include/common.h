#ifndef common_h
#define common_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TTip.h"
#include "TParserLibrary.h"
#include "TEnv.h"

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
static Double_t tigtigTGate[2] = {-60, 60}; // TIGRESS - TIGRESS timing window (ns)
static Double_t tiptipTGate[2] = {-200, 200}; // TIP - TIP fit timing window (ns)
//static Double_t tiptipTGate[2] = {-60, 60}; // narrow TIP - TIP fit timing window (ns)
static Double_t tiptigTGate[2] = {-1480, -1000}; // TIP - TIGRESS timing window (ns)
//static Double_t tiptigTGate[2] = {-1880, -600}; // wide TIP - TIGRESS timing window (ns)
static Double_t tigBGOTGate[2] = {0, 380}; // TIGRESS - BGO timing window (ns)

//PID gates
class PIDGates{
	public :

		PIDGates(); //see common.cxx
		TCutG *alphaRingCut[NTIPRING];
		TCutG *protonRingCut[NTIPRING];
};

static Int_t tip_waveform_pretrigger = 250;

//static Double_t betaCompound = 0.02953; //compound nucleus recoil beta (26Mg)
//static Double_t betaCompound = 0.04143; //compound nucleus recoil beta (31Si)
//static Double_t betaCompound = 0.04093; //compound nucleus recoil beta (32Si)
static Double_t betaCompound = 0.04243; //compound nucleus recoil beta (32Si, high E)
//static Double_t betaCompound = 0.04343; //compound nucleus recoil beta (32Si, 2362 keV line)
static Int_t compoundM_AMU = 33.96786701; //compound mass in atomic mass units (34S)

static Int_t noPileupKValue = 700; //should be 0 for TIG-10s, 700 for GRIF-16s


//FUNCTION PROTOTYPES
Int_t getParticleTypePID(double_t tipPID, double_t energy, Int_t detNum, PIDGates *gates);
Int_t getParticleType(TTipHit *tip_hit, PIDGates *gates);
TVector3 getTigVector(uint8_t core, uint8_t seg);
double_t getEDoppFusEvapManual(double eTig, uint8_t core, uint8_t seg, uint8_t numCsIHits, csi_hit *tip_hits, PIDGates *gates);
double_t getEDoppFusEvapDirect(tig_hit *add_hit, uint8_t numCsIHits, csi_hit *tip_hits, PIDGates *gates);
double_t getEDoppFusEvap(TTigressHit *add_hit, TTip *tip, const uint64_t passedtimeGate, PIDGates *gates);
double_t getTipFitTime(TTipHit *tip_hit, const Int_t pretrigger_samples);
bool ExptSuppression(TDetectorHit* tig, TBgoHit& bgo);
bool gate1D(const Double_t value, const Double_t min, const Double_t max);
Int_t getTIPRing(const Int_t tipPosition);
Int_t getTIGRESSRing(const float theta);
Int_t getTIGRESSSegmentRing(const float theta);
uint64_t passesTimeGateNoAB(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit);
uint64_t passesTimeGateAB(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit);
uint64_t passesTimeGate(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit);

#endif
