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
TH1F *hpgeE, *addE, *addE_ring[NTIGRING], *hpgeRate;
TH2F *hpgeE_ANum, *addE_ANum, *addE_theta, *addE_phi, *theta_phi, *hpgeNum_time;

//Timing
TH1F *hpgeT_hpgeT, *addT_addT;
TH1F *hpgeT_hpgeT_tsep, *addT_addT_tsep;
TH1F *hpgeT_hpgeT_tsepmult2, *addT_addT_tsepmult2;

//HPGe-HPGe
TH2I *hpgePos_hpgePos;
TH2F *hpgeE_hpgeE, *addE_addE;
TH2F *hpgeE_hpgeE_tsep, *addE_addE_tsep;
TH2F *hpgeE_hpgeE_tsepmult2, *addE_addE_tsepmult2;


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
  hpgeNum_time = new TH2F("Tigress_Array_Number_vs_Time","Tigress Array Number vs Time;Time (s);Array Number",3600,0,3600,64,0,64);
  hpgeNum_time->GetXaxis()->SetTitle("Run Time (s)");
  hpgeList->Add(hpgeNum_time);
  hpgeE = new TH1F("HPGe_Energy", "Tigress Energy (non-addback)", 16384, 0, 8192);
  hpgeList->Add(hpgeE);
  hpgeE_ANum = new TH2F("HPGe_Energy_vs_Array_Number", "HPGe Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  hpgeList->Add(hpgeE_ANum);
  addE = new TH1F("Addback_Energy", "Addback Energy", 16384, 0, 8192);
  hpgeList->Add(addE);
  addE_ANum = new TH2F("Addback_Energy_vs_Array_Number", "Addback Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  hpgeList->Add(addE_ANum);
  for(int i=0; i<NTIGRING; i++){
    addE_ring[i] = new TH1F(Form("Addback_energy_ring_%i",i+1), Form("Addback energy (ring %i)",i+1), 8192, 0, 8192);
    addE_ring[i]->GetXaxis()->SetTitle("Addback Energy");
    hpgeList->Add(addE_ring[i]);
  }
  addE_theta = new TH2F("Addback_Energy_vs_theta_segment", "Addback Energy vs. theta (segment)", 1800, 0, 180, 4096, 0, 8192);
  addE_theta->GetXaxis()->SetTitle("#theta (deg)");
  addE_theta->GetYaxis()->SetTitle("Addback Energy");
  hpgeList->Add(addE_theta);
  addE_phi = new TH2F("Addback_Energy_vs_phi_segment", "Addback Energy vs. phi (segment)", 1800, -180, 180, 4096, 0, 8192);
  addE_phi->GetXaxis()->SetTitle("#phi (deg)");
  addE_phi->GetYaxis()->SetTitle("Addback Energy");
  hpgeList->Add(addE_phi);
  theta_phi = new TH2F("theta_vs_phi_segment", "theta vs. phi (segment, addback data)", 1800, 0, 180, 1800, -180, 180);
  theta_phi->GetXaxis()->SetTitle("#theta (deg)");
  theta_phi->GetYaxis()->SetTitle("#phi (deg)");
  hpgeList->Add(theta_phi);
  hpgeRate = new TH1F("HPGe_Total_Rate", "HPGe Total Rate", 8192, 0, 8192);
  hpgeRate->GetXaxis()->SetTitle("Run Time (s)");
  hpgeRate->GetYaxis()->SetTitle("Count/s");
  hpgeList->Add(hpgeRate);

  //Timing spectra
  hpgeT_hpgeT = new TH1F("Tigress_Tigress_time","Tigress - Tigress time",4096,-2048,2048); 
  hpgeT_hpgeT->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT);
  addT_addT = new TH1F("Addback_Addback_time","Addback - Addback time",4096,-2048,2048); 
  addT_addT->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(addT_addT);
  hpgeT_hpgeT_tsep = new TH1F("Tigress_Tigress_time_tsep","Tigress - Tigress time, time separated",4096,-2048,2048); 
  hpgeT_hpgeT_tsep->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsep);
  addT_addT_tsep = new TH1F("Addback_Addback_time_tsep","Addback - Addback time, time separated",4096,-2048,2048); 
  addT_addT_tsep->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(addT_addT_tsep);
  hpgeT_hpgeT_tsepmult2 = new TH1F("Tigress_Tigress_time_tsepmult2","Tigress - Tigress time, time separated multiplicity 2",4096,-2048,2048); 
  hpgeT_hpgeT_tsepmult2->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsepmult2);
  addT_addT_tsepmult2 = new TH1F("Addback_Addback_time_tsepmult2","Addback - Addback time, time separated multiplicity 2",4096,-2048,2048); 
  addT_addT_tsepmult2->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(addT_addT_tsepmult2);

  //HPGe-HPGe
  hpgePos_hpgePos = new TH2I("HPGe_HPGe_position", "HPGe-HPGe position", 64, 0, 64, 64, 0, 64);
  hpgePos_hpgePos->GetXaxis()->SetTitle("HPGe position 1");
  hpgePos_hpgePos->GetYaxis()->SetTitle("HPGe position 2");
  hpgehpgeList->Add(hpgePos_hpgePos);
  hpgeE_hpgeE = new TH2F("HPGe_Gamma_Gamma", "HPGe Gamma-Gamma", 4096, 0, 8192, 4096, 0, 8192);
  hpgeE_hpgeE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE);
  addE_addE = new TH2F("Addback_Gamma_Gamma", "Addback Gamma-Gamma", 4096, 0, 8192, 4096, 0, 8192);
  addE_addE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(addE_addE);
  hpgeE_hpgeE_tsep = new TH2F("HPGe_Gamma_Gamma_tsep", "HPGe Gamma-Gamma time separated", 4096, 0, 8192, 4096, 0, 8192);
  hpgeE_hpgeE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsep);
  addE_addE_tsep = new TH2F("Addback_Gamma_Gamma_tsep", "Addback Gamma-Gamma time separated", 4096, 0, 8192, 4096, 0, 8192);
  addE_addE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(addE_addE_tsep);
  hpgeE_hpgeE_tsepmult2 = new TH2F("HPGe_Gamma_Gamma_tsepmult2", "HPGe Gamma-Gamma time separated multiplicity 2", 4096, 0, 8192, 4096, 0, 8192);
  hpgeE_hpgeE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsepmult2);
  addE_addE_tsepmult2 = new TH2F("Addback_Gamma_Gamma_tsepmult2", "Addback Gamma-Gamma time separated multiplicity 2", 4096, 0, 8192, 4096, 0, 8192);
  addE_addE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(addE_addE_tsepmult2);

}
