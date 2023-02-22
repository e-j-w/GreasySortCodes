#ifndef SortDiagnosticsG_h
#define SortDiagnosticsG_h

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TTigress.h"
#include "TTip.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

TApplication *theApp;

TList *tigList, *tipList, *timingList, *tipPIDList, *tipPIDGateList, *tiptipList, *tigtigList;
TList *tiptigList, *tigPIDSepList, *tigtigPIDSepList; 

//Raw TIGRESS
TH1F *addE, *addE_ring[NTIGRING], *tigRate;
TH2F *addE_ANum, *addE_theta, *tigNum_time;

//TIP
TH1F *tip_E, *tip_Etot, *tipRate;
TH2F *tip_E_pos;
TH1I *tip_mult, *tip_pos, *tip_ring, *tip_fittype;

//Timing
TH1F *tiptipFitT;
TH1F *addT_addT;
TH1F *tipT_tigT_diff;
TH2F *tipTtigT_addE, *tipTtigT_EDopp, *tipTtigT_EDopp_no90, *tigTtigT_addE;

//TIP PID
TH2F *tip_E_PID_Sum, *tip_E_PID_Ring[NTIPRING], *tip_E_PID[NTIP];

//TIP-TIP
TH2I *tipPos_tipPos;

//TIGRESS-TIGRESS
TH2I *tigPos_tigPos;
TH2F *addE_addE;

//TIGRESS-TIP
TH2I *tiptig_Amult;
TH2F *tigE_tipTtigTdiff, *tipE_tipTStigTdiff, *tipPos_tipTStigTdiff, *tigE_tipRing, *tigE_tipMult, *tigE_tipE;

//TIGRESS PID Separated plots
TH2F *addE_xayp_ring[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH2F *addDopp_xayp_ring[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

//TIGRESS-TIGRESS PID Separated plots
TH2F *addEaddE_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH2F *addDoppaddDopp_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];
TH2F *addEaddDopp_xayp[MAX_NUM_PARTICLE+1][MAX_NUM_PARTICLE+1];

PIDGates *gates;

class SortDiagnostics {

	public :

		SortDiagnostics(){;} 
		void SortData(const char*, const char*);
		void Initialise();
};
#endif

void SortDiagnostics::Initialise() {

  cout << "Creating lists" << endl;

  tigList = new TList;
  tipList = new TList;
  timingList = new TList;
  tipPIDList = new TList;
  tipPIDGateList = new TList;
  tiptipList = new TList;
  tigtigList = new TList;
  tiptigList = new TList;
  tigPIDSepList = new TList;
  tigtigPIDSepList = new TList;

  cout << "Creating histograms" << endl;
  //histogram names shouldn't use spaces, to aid GRSISort-based analysis

  //Raw TIGRESS Spectra
  tigNum_time = new TH2F("Tigress_Array_Number_vs_Time","Tigress Array Number vs Time;Time (s);Array Number",3600,0,3600,64,0,64);
  tigNum_time->GetXaxis()->SetTitle("Run Time (s)");
  tigList->Add(tigNum_time);
  addE = new TH1F("Addback_Energy", "Addback Energy", 16384, 0, 8192);
  tigList->Add(addE);
  addE_ANum = new TH2F("Addback_Energy_vs_Array_Number", "Addback Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(addE_ANum);
  for(int i=0; i<NTIGRING; i++){
    addE_ring[i] = new TH1F(Form("Addback_energy_ring_%i",i+1), Form("Addback energy (ring %i)",i+1), 8192, 0, 8192);
    addE_ring[i]->GetXaxis()->SetTitle("Addback Energy");
    tigList->Add(addE_ring[i]);
  }
  addE_theta = new TH2F("Addback_Energy_vs_theta_segment", "Addback Energy vs. theta (segment)", 1800, 0, 180, 8192, 0, 8192);
  addE_theta->GetXaxis()->SetTitle("#theta (deg)");
  addE_theta->GetYaxis()->SetTitle("Addback Energy");
  tigList->Add(addE_theta);
  tigRate = new TH1F("TIGRESS_Total_Rate", "TIGRESS Total Rate", 8192, 0, 8192);
  tigRate->GetXaxis()->SetTitle("Run Time (s)");
  tigRate->GetYaxis()->SetTitle("Count/s");
  tigList->Add(tigRate);

  //Raw TIP Spectra
  tip_pos = new TH1I("TIP_position","TIP position",129,0,129); 
  tipList->Add(tip_pos);
  tip_ring = new TH1I("TIP_ring","TIP ring",10,0,10);
  tipList->Add(tip_ring);
  tip_E = new TH1F("TIP_energy","TIP energy",8192,0,128); 
  tipList->Add(tip_E);
  tip_Etot = new TH1F("TIP_energy_total_per_event","TIP energy (total per-event)",8192,0,128); 
  tipList->Add(tip_Etot);
  tip_mult = new TH1I("TIP_multiplicity","TIP multiplicity",10,0,10); 
  tipList->Add(tip_mult);
  tip_fittype = new TH1I("TIP waveform fit type","TIP waveform fit type",5,0,5);
  tipList->Add(tip_fittype);
  tip_E_pos = new TH2F("TIP_energy_vs_position","TIP energy vs position",8192,0,128,129,0,129);
  tip_E_pos->GetXaxis()->SetTitle("E (arb.)");
  tip_E_pos->GetYaxis()->SetTitle("TIP position");
  tipList->Add(tip_E_pos);
  tipRate = new TH1F("TIP_Total_Rate", "TIP Total Rate", 8192, 0, 8192);
  tipRate->GetXaxis()->SetTitle("Run Time (s)");
  tipRate->GetYaxis()->SetTitle("Count/s");
  tipList->Add(tipRate);

  //Timing spectra
  tiptipFitT = new TH1F("TIP_TIP_time_fit","TIP - TIP time (fit)",8192,-4096,4096);
  tiptipFitT->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIP, fit} (ns)");
  timingList->Add(tiptipFitT);
  addT_addT = new TH1F("Tigress_Tigress_time","Tigress - Tigress time",4096,-2048,2048); 
  addT_addT->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(addT_addT);
  tipT_tigT_diff = new TH1F("TIP_fit_Tigress_time","TIP fit - Tigress time",4096,-4096,4096);
  tipT_tigT_diff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  timingList->Add(tipT_tigT_diff);
  tipTtigT_addE = new TH2F("addE_TIP_fit_Tigress_time","Addback energy vs. Tigress - TIP (avg of hits) fit time",4096,0,8192,4096,0,4096);
  tipTtigT_addE->GetYaxis()->SetTitle("t_{TIGRESS} (ns) - t_{avg TIP, fit}");
  tipTtigT_addE->GetXaxis()->SetTitle("Addback Energy");
  timingList->Add(tipTtigT_addE);
  tipTtigT_EDopp = new TH2F("EDopp_TIP_fit_Tigress_time","Doppler energy vs. Tigress - TIP (avg of hits) fit time",4096,0,8192,4096,0,4096);
  tipTtigT_EDopp->GetYaxis()->SetTitle("t_{TIGRESS} (ns) - t_{avg TIP, fit}");
  tipTtigT_EDopp->GetXaxis()->SetTitle("Doppler Energy");
  timingList->Add(tipTtigT_EDopp);
  tipTtigT_EDopp_no90 = new TH2F("EDopp_no90_TIP_fit_Tigress_time","Doppler energy (no 90 degree positions) vs. Tigress - TIP (avg of hits) fit time",4096,0,8192,4096,0,4096);
  tipTtigT_EDopp_no90->GetYaxis()->SetTitle("t_{TIGRESS} - t_{avg TIP, fit} (ns)");
  tipTtigT_EDopp_no90->GetXaxis()->SetTitle("Doppler Energy");
  timingList->Add(tipTtigT_EDopp_no90);
  tigTtigT_addE = new TH2F("addE_Tigress-Tigress_time","Doppler energy vs. Tigress - Tigress time",4096,0,8192,4096,-4096,4096);
  tigTtigT_addE->GetYaxis()->SetTitle("t_{TIGRESS,2} - t_{TIGRESS,1} (ns)");
  tigTtigT_addE->GetXaxis()->SetTitle("Addback Energy");
  timingList->Add(tigTtigT_addE);

  //TIP PID
  tip_E_PID_Sum = new TH2F("TIP_energy_vs_PID_Sum","TIP energy vs PID (Sum)",1024,0,128,512,0,512);
  tip_E_PID_Sum->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
  tip_E_PID_Sum->GetXaxis()->SetTitle("CsI E (MeV)");
  tipPIDList->Add(tip_E_PID_Sum);
  for(int i=0; i<NTIPRING; i++){
    tip_E_PID_Ring[i] = new TH2F(Form("TIP_energy_vs_PID_Ring_%i",i),Form("TIP energy vs PID (Ring %i)",i),1024,0,128,512,0,512);
    tip_E_PID_Ring[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
    tip_E_PID_Ring[i]->GetXaxis()->SetTitle("CsI E (MeV)");
    tipPIDList->Add(tip_E_PID_Ring[i]);
  }
  for(int i=0; i<NTIP; i++){
    tip_E_PID[i] = new TH2F(Form("TIP_energy_vs_PID_Pos_%i",i+1),Form("TIP energy vs PID (Pos %i)",i+1),1024,0,128,512,0,512);
    tip_E_PID[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100");
    tip_E_PID[i]->GetXaxis()->SetTitle("CsI E (MeV)");
    tipPIDList->Add(tip_E_PID[i]);
  }

  //TIP-TIP
  tipPos_tipPos = new TH2I("TIP_TIP_position", "TIP-TIP position", 129, 0, 129, 129, 0, 129);
  tipPos_tipPos->GetXaxis()->SetTitle("TIP position 1");
  tipPos_tipPos->GetYaxis()->SetTitle("TIP position 2");
  tiptipList->Add(tipPos_tipPos);

  //TIGRESS-TIGRESS
  tigPos_tigPos = new TH2I("TIGRESS_TIGRESS_position", "TIGRESS-TIGRESS position", 64, 0, 64, 64, 0, 64);
  tigPos_tigPos->GetXaxis()->SetTitle("TIGRESS position 1");
  tigPos_tigPos->GetYaxis()->SetTitle("TIGRESS position 2");
  tigtigList->Add(tigPos_tigPos);
  addE_addE = new TH2F("Addback_Gamma_Gamma", "Addback Gamma-Gamma", 4096, 0, 8192, 4096, 0, 8192);
  addE_addE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(addE_addE);

  //TIGRESS-TIP
  tiptig_Amult = new TH2I("TIP_TIGRESS_addback_multiplicity","TIP-TIGRESS addback multiplicity",10,0,10,10,0,10);
  tiptig_Amult->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_Amult->GetYaxis()->SetTitle("TIGRESS multiplicity (unsuppressed addback)");
  tiptigList->Add(tiptig_Amult);
  tigE_tipTtigTdiff = new TH2F("Tigress_addback_energy_vs_TIP_fit_Tigress_time", "Tigress addback energy vs TIP fit-Tigress time",4096,-4096,4096,8192,0,8192);
  tigE_tipTtigTdiff->GetYaxis()->SetTitle("TIGRESS addback energy (keV)");
  tigE_tipTtigTdiff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tigE_tipTtigTdiff);
  tigE_tipRing = new TH2F("Tigress_addback_energy_vs_TIP_ring", "Tigress addback energy vs TIP ring", 8192, 0, 8192, 10, 0, 10);
  tigE_tipRing->GetYaxis()->SetTitle("TIP ring");
  tigE_tipRing->GetXaxis()->SetTitle("TIGRESS energy");
  tiptigList->Add(tigE_tipRing);
  tigE_tipMult = new TH2F("Tigress_addback_energy_vs_TIP_multiplicity", "Tigress addback energy vs TIP multiplicity", 8192, 0, 8192, 10, 0, 10);
  tigE_tipMult->GetYaxis()->SetTitle("TIP multiplicity");
  tigE_tipMult->GetXaxis()->SetTitle("TIGRESS energy");
  tiptigList->Add(tigE_tipMult);
  tigE_tipE = new TH2F("Tigress_addback_energy_vs_TIP_energy", "Tigress addback energy vs TIP energy;TIP Energy;Tigress addback Energy",2000,0,4000,4000,0,4000);
  tiptigList->Add(tigE_tipE);
  tipE_tipTStigTdiff = new TH2F("TIP_TS_Tigress_time_vs_TIP_energy","TIP TS-Tigress time vs TIP energy",4096,-4096,4096,2048,0,8192);
  tipE_tipTStigTdiff->GetYaxis()->SetTitle("TIP E (arb.)");
  tipE_tipTStigTdiff->GetXaxis()->SetTitle("t_{TIP, TS} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipE_tipTStigTdiff);
  tipPos_tipTStigTdiff = new TH2F("TIP_TS_Tigress_time_vs_TIP_position","TIP TS-Tigress time vs TIP position",4096,-4096,4096,129,0,129);
  tipPos_tipTStigTdiff->GetYaxis()->SetTitle("TIP position");
  tipPos_tipTStigTdiff->GetXaxis()->SetTitle("t_{TIP, TS} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipPos_tipTStigTdiff);

  //TIGRESS PID separated
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS ring spectra
        addE_xayp_ring[i][j] = new TH2F(Form("TIGRESS_addback_energy_%ip%ia_gate",i,j),Form("TIGRESS addback energy (%ip%ia gate)",i,j),8192,0,8192,NTIGSEGRING+1,0,NTIGSEGRING+1);
        addE_xayp_ring[i][j]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        addE_xayp_ring[i][j]->GetYaxis()->SetTitle("TIGRESS Ring (ring 0 = sum)");
        tigPIDSepList->Add(addE_xayp_ring[i][j]);
      }
    }
  }
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS ring spectra
        addDopp_xayp_ring[i][j] = new TH2F(Form("Doppler_corrected_energy_%ip%ia_gate",i,j),Form("TIGRESS Doppler corrected addback energy (%ip%ia gate, beta=%f)",i,j,betaCompound),8192,0,8192,NTIGSEGRING+1,0,NTIGSEGRING+1);
        addDopp_xayp_ring[i][j]->GetXaxis()->SetTitle("E_{#gamma} (keV)");
        addDopp_xayp_ring[i][j]->GetYaxis()->SetTitle("TIGRESS Ring (ring 0 = sum)");
        tigPIDSepList->Add(addDopp_xayp_ring[i][j]);
      }
    }
  }

  //TIGRESS-TIGRESS PID separated
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS ring spectra
        addEaddE_xayp[i][j] = new TH2F(Form("Addback_Gamma_Gamma_%ip%ia_gate",i,j),Form("Addback Gamma-Gamma (%ip%ia gate)",i,j),8192,0,8192,8192,0,8192);
        addEaddE_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
        addEaddE_xayp[i][j]->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
        tigtigPIDSepList->Add(addEaddE_xayp[i][j]);
      }
    }
  }
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS ring spectra
        addDoppaddDopp_xayp[i][j] = new TH2F(Form("Doppler_corrected_Gamma_Gamma_%ip%ia_gate",i,j),Form("Doppler corrected Gamma-Gamma (%ip%ia gate, beta=%f)",i,j,betaCompound),8192,0,8192,8192,0,8192);
        addDoppaddDopp_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma 1} (keV, Doppler corrected)");
        addDoppaddDopp_xayp[i][j]->GetYaxis()->SetTitle("E_{#gamma 2} (keV, Doppler corrected)");
        tigtigPIDSepList->Add(addDoppaddDopp_xayp[i][j]);
      }
    }
  }
  for(int i=0; i<MAX_NUM_PARTICLE+1; i++){
    for(int j=0; j<MAX_NUM_PARTICLE+1; j++){
      if((i+j)<=MAX_NUM_PARTICLE){
        //TIGRESS ring spectra
        addEaddDopp_xayp[i][j] = new TH2F(Form("Addback_Gamma_Doppler_corrected_Gamma_%ip%ia_gate",i,j),Form("Addback Gamma-Doppler corrected Gamma (%ip%ia gate, beta=%f)",i,j,betaCompound),8192,0,8192,8192,0,8192);
        addEaddDopp_xayp[i][j]->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
        addEaddDopp_xayp[i][j]->GetYaxis()->SetTitle("E_{#gamma 2} (keV, Doppler corrected)");
        tigtigPIDSepList->Add(addEaddDopp_xayp[i][j]);
      }
    }
  }

  //Setup TIP PID gates
  gates = new PIDGates;
  for(int i=0; i<NTIPRING; i++){
    tipPIDGateList->Add(gates->alphaRingCut[i]);
    tipPIDGateList->Add(gates->protonRingCut[i]);
  }

}
