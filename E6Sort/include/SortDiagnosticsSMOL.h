#ifndef SortDiagnosticsS_h
#define SortDiagnosticsS_h

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"

#include "TCanvas.h"
#include "TApplication.h"

#include "evt_fmt.h"
#include <stdint.h> //allows uint8_t and similiar types

using namespace std;

TApplication *theApp;

TList *hpgeList, *tipList, *timingList, *tipPIDList, *tipPIDGateList, *tiptipList, *hpgehpgeList;
TList *tiphpgeList, *hpgePIDSepList, *hpgehpgePIDSepList; 

//Raw HPGe
TH1F *hpgeE, *hpgeRate;
TH2F *hpgeE_ANum, *hpgeNum_time;
TH1I *hpgeMult;

//Timing
TH1F *hpgeT_hpgeT;
TH1F *hpgeT_hpgeT_tsep;
TH1F *hpgeT_hpgeT_tsepmult2;
TH1F *hpgeT_hpgeT_tseprand;

//HPGe-HPGe
TH1F *hpge_hpge_dist, *hpge_hpge_angle;
TH2I *hpgePos_hpgePos;
TH2F *hpgeE_hpgeE;
TH2F *hpgeE_hpgeE_tsep;
TH2F *hpgeE_hpgeE_tsepmult2;
TH2F *hpgeE_hpgeE_tseprand;
TH2F *hpgeE_hpgeE_180deg;
TH1F *hpgeE_hpgeE_180deg_proj;
TH1F *hpgeE_hpgeE_180deg_sum;

class SortDiagnosticsS {

	public :

		SortDiagnosticsS(){;} 
		void SortData(const char*, const char*);
		void Initialise();
};
#endif

void SortDiagnosticsS::Initialise() {

  cout << "Creating lists" << endl;

  hpgeList = new TList;
  timingList = new TList;
  hpgehpgeList = new TList;
  hpgePIDSepList = new TList;
  hpgehpgePIDSepList = new TList;

  cout << "Creating histograms" << endl;
  //histogram names shouldn't use spaces, to aid GRSISort-based analysis

  //Raw HPGe Spectra
  hpgeNum_time = new TH2F("HPGe_Array_Number_vs_Time","HPGe Array Number vs Time;Time (s);Array Number",3600,0,3600,64,0,64);
  hpgeNum_time->GetXaxis()->SetTitle("Run Time (s)");
  hpgeList->Add(hpgeNum_time);
  hpgeE = new TH1F("HPGe_Energy", "HPGe Energy (non-addback)", 16384, 0, 8192);
  hpgeList->Add(hpgeE);
  hpgeE_ANum = new TH2F("HPGe_Energy_vs_Array_Number", "HPGe Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  hpgeList->Add(hpgeE_ANum);
  hpgeRate = new TH1F("HPGe_Total_Rate", "HPGe Total Rate", 8192, 0, 8192);
  hpgeRate->GetXaxis()->SetTitle("Run Time (s)");
  hpgeRate->GetYaxis()->SetTitle("Count/s");
  hpgeList->Add(hpgeRate);
  hpgeMult = new TH1I("hpgeMult", "HPGe multiplicity (non-addback)", 64, 0, 64);
  hpgeMult->GetYaxis()->SetTitle("Counts");
  hpgeList->Add(hpgeMult);

  //Timing spectra
  hpgeT_hpgeT = new TH1F("HPGe_HPGe_time","HPGe - HPGe time",4096,-2048,2048); 
  hpgeT_hpgeT->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT);
  hpgeT_hpgeT_tsep = new TH1F("HPGe_HPGe_time_tsep","HPGe - HPGe time, time separated",4096,-2048,2048); 
  hpgeT_hpgeT_tsep->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsep);
  hpgeT_hpgeT_tsepmult2 = new TH1F("HPGe_HPGe_time_tsepmult2","HPGe - HPGe time, time separated multiplicity 2",4096,-2048,2048); 
  hpgeT_hpgeT_tsepmult2->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsepmult2);
  hpgeT_hpgeT_tseprand = new TH1F("HPGe_HPGe_time_tseprand","HPGe - HPGe time, time-random separated",4096,-2048,2048); 
  hpgeT_hpgeT_tseprand->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tseprand);

  //HPGe-HPGe
  hpge_hpge_dist = new TH1F("HPGe_HPGe_distance","HPGe-HPGe distance",512,0,512);
  hpge_hpge_dist->GetXaxis()->SetTitle("HPGe - HPGe distance (mm)");
  hpgehpgeList->Add(hpge_hpge_dist);
  hpge_hpge_angle = new TH1F("HPGe_HPGe_angle","HPGe-HPGe angle",512,0,360);
  hpge_hpge_angle->GetXaxis()->SetTitle("HPGe - HPGe angle (deg)");
  hpgehpgeList->Add(hpge_hpge_angle);
  hpgePos_hpgePos = new TH2I("HPGe_HPGe_position", "HPGe-HPGe position", 64, 0, 64, 64, 0, 64);
  hpgePos_hpgePos->GetXaxis()->SetTitle("HPGe position 1");
  hpgePos_hpgePos->GetYaxis()->SetTitle("HPGe position 2");
  hpgehpgeList->Add(hpgePos_hpgePos);
  hpgeE_hpgeE = new TH2F("HPGe_Gamma_Gamma", "HPGe Gamma-Gamma", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE);
  hpgeE_hpgeE_tsep = new TH2F("HPGe_Gamma_Gamma_tsep", "HPGe Gamma-Gamma time separated", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsep);
  hpgeE_hpgeE_tsepmult2 = new TH2F("HPGe_Gamma_Gamma_tsepmult2", "HPGe Gamma-Gamma time separated multiplicity 2", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsepmult2);
  hpgeE_hpgeE_tseprand = new TH2F("HPGe_Gamma_Gamma_tseprand", "HPGe Gamma-Gamma time-random separated", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tseprand->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tseprand->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tseprand);
  hpgeE_hpgeE_180deg = new TH2F("HPGe_Gamma_Gamma_180deg", "HPGe Gamma-Gamma 180 degree coincidences", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_180deg->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_180deg->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg);
  hpgeE_hpgeE_180deg_proj = new TH1F("HPGe_Gamma_Gamma_180deg_proj", "HPGe Gamma-Gamma 180 degree coincidences (projection)", 8192, 0, 4096);
  hpgeE_hpgeE_180deg_proj->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_proj);
  hpgeE_hpgeE_180deg_sum = new TH1F("HPGe_Gamma_Gamma_180deg_sum", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits)", 8192, 0, 4096);
  hpgeE_hpgeE_180deg_sum->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum);

}
