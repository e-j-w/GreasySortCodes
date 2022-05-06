#ifndef SortCode_h
#define SortCode_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCutG.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TReaction.h"
#include "TSRIM.h"
#include "TTigress.h"
#include "TTip.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TEmma.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;

TList *tigList, *tipList, *tipPIDList, *tipPIDGateList, *tigtigList, *tigbgoList, *tiptigList; 

//Raw TIGRESS
TH1F *tigE, *tigE_unsupp, *addE, *addE_ring[NTIGRING];
TH2F *tigE_ANum,*addE_ANum, *addE_theta;
TH1F *tigChan;

//TIP
TH1F *tip_E, *tip_Etot, *tip_CFDFitDiff, *tip_wfrmsize, *tiptipT, *tiptipFitT;
TH2F *tip_E_pos, *tip_E_waveformAmp;
TH1I *tip_mult, *tip_pos, *tip_ring;

//TIP PID
TH2F *tip_E_PID_Sum, *tip_E_PID_Ring[NTIPRING], *tip_E_PID[NTIP];

//TIGRESS-TIGRESS
TH1F *addT_addT;
TH2F *addE_addE, *addE_addE_tg;

//TIGRESS-BGO
TH1I *bgo_mult, *bgo_det;
TH1F *tigT_bgoT, *tigT_bgoT_supp;

//TIGRESS-TIP
TH2I *tiptig_mult, *tiptig_multSupp;
TH1F *tipT_tigT_diff, *tipTCFD_tigT_diff, *tigE_TIPtg;
TH2F *tigE_tipTtigTdiff, *tigE_tipRing;
TH2F *tigE_tipMult;

//PID gates
TCutG *alphaRingCut[NTIPRING], *protonRingCut[NTIPRING];

class SortCode {

	public :

		SortCode(){;} 
		void SortData(const char*, const char*, const char*);
		void Initialise();
};
#endif

void SortCode::Initialise() {

  printf("Start initialization\n");
  printf("Creating lists\n");

  tigList = new TList;
  tipList = new TList;
  tipPIDList = new TList;
  tipPIDGateList = new TList;
  tigtigList = new TList;
  tigbgoList = new TList;
  tiptigList = new TList;

  printf("Creating histograms\n");

  //Raw TIGRESS Spectra
  tigE_unsupp = new TH1F("Tigress Energy (unsuppressed)", "Tigress Energy (unsuppressed)", 8192, 0, 8192);
  tigList->Add(tigE_unsupp);
  tigE = new TH1F("Tigress Energy", "Tigress Energy", 8192, 0, 8192);
  tigList->Add(tigE);
  tigE_ANum = new TH2F("Tigress Energy vs. Array Number", "Tigress Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigE_ANum);
  addE = new TH1F("Addback Energy", "Addback Energy", 8192, 0, 8192);
  tigList->Add(addE);
  addE_ANum = new TH2F("Addback Energy vs. Array Number", "Addback Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(addE_ANum);
  for(int i=0; i<NTIGRING; i++){
    addE_ring[i] = new TH1F(Form("Addback energy (ring %i)",i+1), Form("Addback energy (ring %i)",i+1), 8192, 0, 8192);
    addE_ring[i]->GetXaxis()->SetTitle("Addback Energy");
    tigList->Add(addE_ring[i]);
  }
  addE_theta = new TH2F("Addback Energy vs. theta (segment)", "Addback Energy vs. theta (segment)", 180, 0, 180, 8192, 0, 8192);
  addE_theta->GetXaxis()->SetTitle("#{theta} (deg)");
  addE_theta->GetYaxis()->SetTitle("Addback Energy");
  tigList->Add(addE_theta);
  tigChan = new TH1F("TIGRESS Channel Number","Tigress Channel Number;Channel;Counts/channel",64,0,64);
  tigList->Add(tigChan);

  //Raw TIP Spectra
  tip_pos = new TH1I("TIP position","TIP position",129,0,129); 
  tipList->Add(tip_pos);
  tip_ring = new TH1I("TIP ring","TIP ring",10,0,10);
  tipList->Add(tip_ring);
  tip_E = new TH1F("TIP energy","TIP energy",8192,0,128); 
  tipList->Add(tip_E);
  tip_Etot = new TH1F("TIP energy (total per-event)","TIP energy (total per-event)",8192,0,128); 
  tipList->Add(tip_Etot);
  tip_mult = new TH1I("TIP multiplicity","TIP multiplicity",10,0,10); 
  tipList->Add(tip_mult);
  tip_wfrmsize = new TH1F("TIP waveform size","TIP waveform size",4096,0,4096); 
  tipList->Add(tip_wfrmsize);
  tip_CFDFitDiff = new TH1F("TIP Fit - CFD timing","TIP Fit - CFD timing",2048,-1024,1024);
  tip_CFDFitDiff->GetXaxis()->SetTitle("t_{fit} - t_{CFD} (ns)");
  tipList->Add(tip_CFDFitDiff);
  tip_E_pos = new TH2F("TIP energy vs position","TIP energy vs position",8192,0,128,129,0,129);
  tip_E_pos->GetXaxis()->SetTitle("E (arb.)");
  tip_E_pos->GetYaxis()->SetTitle("TIP position");
  tipList->Add(tip_E_pos);
  tip_E_waveformAmp = new TH2F("TIP energy vs wavform amplitude","TIP energy vs wavform amplitude",4096,0,128,4096,0,8192);
  tip_E_waveformAmp->GetXaxis()->SetTitle("E (arb.)");
  tip_E_waveformAmp->GetYaxis()->SetTitle("Waveform Amplitude (arb.)");
  tipList->Add(tip_E_waveformAmp);
  tiptipT = new TH1F("TIP-TIP timing (CFD)","TIP-TIP timing (CFD)",8192,-4096,4096);
  tipList->Add(tiptipT);
  tiptipFitT = new TH1F("TIP-TIP timing (fit)","TIP-TIP timing (fit)",8192,-4096,4096);
  tipList->Add(tiptipFitT);

  //TIP PID
  tip_E_PID_Sum = new TH2F("TIP energy vs PID (Sum)","TIP energy vs PID (Sum)",1024,0,128,512,0,512);
  tip_E_PID_Sum->GetYaxis()->SetTitle("A_{S}/A_{F} x 100 + 100");
  tip_E_PID_Sum->GetXaxis()->SetTitle("Alpha E (MeV)");
  tipPIDList->Add(tip_E_PID_Sum);
  for(int i=0; i<NTIPRING; i++){
    tip_E_PID_Ring[i] = new TH2F(Form("TIP energy vs PID (Ring %i)",i+1),Form("TIP energy vs PID (Ring %i)",i+1),1024,0,128,512,0,512);
    tip_E_PID_Ring[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100 + 100");
    tip_E_PID_Ring[i]->GetXaxis()->SetTitle("Alpha E (MeV)");
    tipPIDList->Add(tip_E_PID_Ring[i]);
  }
  for(int i=0; i<NTIP; i++){
    tip_E_PID[i] = new TH2F(Form("TIP energy vs PID (Pos %i)",i+1),Form("TIP energy vs PID (Pos %i)",i+1),1024,0,128,512,0,512);
    tip_E_PID[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100 + 100");
    tip_E_PID[i]->GetXaxis()->SetTitle("Alpha E (MeV)");
    tipPIDList->Add(tip_E_PID[i]);
  }

  //TIGRESS-TIGRESS
  addT_addT = new TH1F("addT_addT","Tigress-Tigress_time",4096,-2048,2048); 
  tigtigList->Add(addT_addT);
  addE_addE = new TH2F("addE_addE", "Addback_Gamma_Gamma", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE);
  addE_addE_tg = new TH2F("addE_addE_tg", "Addback_Gamma_Gamma_Time_Gated", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE_tg);

  //TIGRESS-BGO
  bgo_mult = new TH1I("BGO multiplicity","BGO multiplicity",100,0,100); 
  tigbgoList->Add(bgo_mult);
  bgo_det = new TH1I("BGO detector position","BGO detector position",100,0,100); 
  tigbgoList->Add(bgo_det);
  tigT_bgoT = new TH1F("Tigress-BGO time (unsuppressed)","Tigress-BGO time (unsuppressed)",4096,-2048,2048);
  tigT_bgoT->GetXaxis()->SetTitle("Same-detector Tigress-BGO time (ns)");
  tigbgoList->Add(tigT_bgoT);
  tigT_bgoT_supp = new TH1F("Tigress-BGO time (suppressed)","Tigress-BGO time (suppressed)",4096,-2048,2048);
  tigT_bgoT_supp->GetXaxis()->SetTitle("Same-detector Tigress-BGO time (ns)");
  tigbgoList->Add(tigT_bgoT_supp);

  //TIGRESS-TIP
  tiptig_mult = new TH2I("TIP-TIGRESS multiplicity","TIP-TIGRESS multiplicity",10,0,10,10,0,10);
  tiptig_mult->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_mult->GetYaxis()->SetTitle("TIGRESS multiplicity (unsuppressed)");
  tiptigList->Add(tiptig_mult);
  tiptig_multSupp = new TH2I("TIP-TIGRESS multiplicity (suppressed)","TIP-TIGRESS multiplicity (suppressed)",10,0,10,10,0,10);
  tiptig_multSupp->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_multSupp->GetYaxis()->SetTitle("TIGRESS multiplicity (suppressed)");
  tiptigList->Add(tiptig_multSupp);
  tipT_tigT_diff = new TH1F("TIP fit - Tigress time","TIP fit - Tigress time",4096,-4096,4096);
  tipT_tigT_diff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tigE_tipTtigTdiff = new TH2F("Tigress addback energy vs TIP(fit)-Tigress time", "Tigress addback energy vs TIP(fit)-Tigress time",4096,-4096,4096,8192,0,8192);
  tigE_tipTtigTdiff->GetYaxis()->SetTitle("TIGRESS addback energy (keV)");
  tigE_tipTtigTdiff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tigE_tipTtigTdiff);
  tipTCFD_tigT_diff = new TH1F("TIP CFD - Tigress time","TIP CFD - Tigress time", 4096,-4096,4096);
  tipTCFD_tigT_diff->GetXaxis()->SetTitle("t_{TIP, CFD} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipTCFD_tigT_diff);

  tigE_TIPtg = new TH1F("Tigress energy (fit time gated)","Tigress energy (CFD time gated)", 8192, 0, 8192); 
  tiptigList->Add(tigE_TIPtg);

  tigE_tipRing = new TH2F("Tigress addback energy vs TIP ring", "Tigress addback energy vs TIP ring", 8192, 0, 8192, 10, 0, 10);
  tigE_tipRing->GetYaxis()->SetTitle("TIP ring");
  tigE_tipRing->GetXaxis()->SetTitle("TIGRESS energy");
  tiptigList->Add(tigE_tipRing);

  tigE_tipMult = new TH2F("Tigress addback energy vs TIP multiplicity", "Tigress addback energy vs TIP multiplicity", 8192, 0, 8192, 10, 0, 10);
  tigE_tipMult->GetYaxis()->SetTitle("TIP multiplicity");
  tigE_tipMult->GetXaxis()->SetTitle("TIGRESS energy");
  tiptigList->Add(tigE_tipMult);

  //Setup PID gates
  printf("Creating PID gates\n");
  //individual gates for each detector
  for(int i=0; i<NTIPRING; i++){
    alphaRingCut[i] = new TCutG(Form("ring %i alpha cut",i+1),7);
    protonRingCut[i] = new TCutG(Form("ring %i proton cut",i+1),7);

    //CHANGE LATER!! this is just proof of concept
    alphaRingCut[i]->SetPoint(0,10,52);
    alphaRingCut[i]->SetPoint(1,10,150);
    alphaRingCut[i]->SetPoint(2,20,197);
    alphaRingCut[i]->SetPoint(3,20,51);
    alphaRingCut[i]->SetPoint(4,18,56);
    alphaRingCut[i]->SetPoint(5,18,27);
    alphaRingCut[i]->SetPoint(6,10,52);

    protonRingCut[i]->SetPoint(0,10,52);
    protonRingCut[i]->SetPoint(1,10,150);
    protonRingCut[i]->SetPoint(2,20,197);
    protonRingCut[i]->SetPoint(3,20,51);
    protonRingCut[i]->SetPoint(4,18,56);
    protonRingCut[i]->SetPoint(5,18,27);
    protonRingCut[i]->SetPoint(6,10,52);
  }

  //add gates to list
  for(int i=0; i<NTIPRING; i++){
    tipPIDGateList->Add(alphaRingCut[i]);
    tipPIDGateList->Add(protonRingCut[i]);
  }

}
