#ifndef common_h
#define common_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TCutG.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TReaction.h"
#include "TTigress.h"
#include "TTip.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

#define NTIP 128 //number of TIP channels
#define NTIPRING 10 //number of TIP rings
#define NTIGRING 6 //number of TIGRESS rings

#define MAX_NUM_PARTICLE 7 //maximum particle multiplicity

#define MAXNUMTIPHIT 16 //maximum number of TIP hits per event
#define MAXNUMTIGHIT 16 //maximum number of TIGRESS hits per event

#define AMU 931.4941 //atomic mass unit, MeV/c^2

#define S32K      32768 //maximum number of channels per spectrum in .mca and .fmca (changing breaks file compatibility)


//GLOBAL VARIABLES
//(static to avoid multiple declaration when linking)

static Double_t pi = TMath::Pi();
static Double_t r2d = TMath::RadToDeg();
static Double_t d2r = TMath::DegToRad();

//timing windows
static Double_t tigtigTGate[2] = {-60, 60}; // TIGRESS - TIGRESS timing window (ns)
static Double_t tiptipTGate[2] = {-200, 200}; // TIP - TIP timing window (ns)
static Double_t tiptigTGate[2] = {2100, 2300}; // TIP - TIGRESS timing window (ns)

//PID gates
class PIDGates{
	public :

		PIDGates(); //see common.cxx
		TCutG *alphaRingCut[NTIPRING];
    TCutG *protonRingCut[NTIPRING];
};

static Int_t tip_waveform_pretrigger = 250;

static Double_t betaCompound = 0.03757; //compound nucleus recoil beta
static Int_t compoundM_AMU = 93.92877; //compound mass in atomic mass units

static Int_t noPileupKValue = 0; //should be 0 for TIG-10s, 700 for GRIF-16s


//FUNCTION PROTOTYPES
Int_t getParticleType(TTipHit *tip_hit, PIDGates *gates);
double_t getEDoppFusEvap(TTigressHit *add_hit, TTip *tip, const uint32_t passedtimeGate, PIDGates *gates);
double_t getTipFitTime(TTipHit *tip_hit, const Int_t pretrigger_samples);
bool S1232Suppression(TDetectorHit* tig, TBgoHit& bgo);
bool gate1D(const Double_t value, const Double_t min, const Double_t max);
Int_t getTIPRing(const Int_t tipPosition);
Int_t getTIGRESSRing(const float theta);
uint32_t passesTimeGate(TTigress *tigress, TTip *tip);

#endif