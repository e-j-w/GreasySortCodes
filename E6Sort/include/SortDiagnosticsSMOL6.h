#ifndef SortDiagnosticsS6_h
#define SortDiagnosticsS6_h

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

TList *tigList, *tipList, *timingList, *tipPIDList, *tipPIDGateList, *tiptipList, *tigtigList;
TList *tiptigList, *tigPIDSepList, *tigtigPIDSepList; 

//Raw TIGRESS
TH1F *tigE, *addE, *tigE_mult2, *addE_ring[NTIGRING], *tigRate;
TH2F *tigE_ANum, *addE_ANum, *addE_theta, *addE_phi, *theta_phi, *tigNum_time;

//Timing
TH1F *tigT_tigT, *addT_addT;
TH1F *tigT_tigT_tsep, *addT_addT_tsep;
TH1F *tigT_tigT_tsepmult2, *addT_addT_tsepmult2;

//TIGRESS-TIGRESS
TH2I *tigPos_tigPos;
TH2F *tigE_tigE, *addE_addE;
TH2F *tigE_tigE_tsep, *addE_addE_tsep;
TH2F *tigE_tigE_tsepmult2, *addE_addE_tsepmult2;
TH2F *tigE_tigE_trand, *addE_addE_trand;


class SortDiagnosticsS6 {

	public :

		SortDiagnosticsS6(){;} 
		void SortData(const char*, const char*);
		void Initialise();
};
#endif

void SortDiagnosticsS6::Initialise() {

  cout << "Creating lists" << endl;

  tigList = new TList;
  timingList = new TList;
  tigtigList = new TList;
  tigPIDSepList = new TList;
  tigtigPIDSepList = new TList;

  cout << "Creating histograms" << endl;
  //histogram names shouldn't use spaces, to aid GRSISort-based analysis

  //Raw TIGRESS Spectra
  tigNum_time = new TH2F("Tigress_Array_Number_vs_Time","Tigress Array Number vs Time;Time (s);Array Number",3600,0,3600,64,0,64);
  tigNum_time->GetXaxis()->SetTitle("Run Time (s)");
  tigList->Add(tigNum_time);
  tigE = new TH1F("TIGRESS_Energy", "Tigress Energy (non-addback)", 16384, 0, 8192);
  tigList->Add(tigE);
  tigE_mult2 = new TH1F("TIGRESS_Energy_mult2", "Tigress Energy (non-addback multiplicity 2)", 16384, 0, 8192);
  tigList->Add(tigE_mult2);
  tigE_ANum = new TH2F("TIGRESS_Energy_vs_Array_Number", "TIGRESS Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigE_ANum);
  addE = new TH1F("Addback_Energy", "Addback Energy", 16384, 0, 8192);
  tigList->Add(addE);
  addE_ANum = new TH2F("Addback_Energy_vs_Array_Number", "Addback Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(addE_ANum);
  for(int i=0; i<NTIGRING; i++){
    addE_ring[i] = new TH1F(Form("Addback_energy_ring_%i",i+1), Form("Addback energy (ring %i)",i+1), 8192, 0, 8192);
    addE_ring[i]->GetXaxis()->SetTitle("Addback Energy");
    tigList->Add(addE_ring[i]);
  }
  addE_theta = new TH2F("Addback_Energy_vs_theta_segment", "Addback Energy vs. theta (segment)", 1800, 0, 180, 4096, 0, 8192);
  addE_theta->GetXaxis()->SetTitle("#theta (deg)");
  addE_theta->GetYaxis()->SetTitle("Addback Energy");
  tigList->Add(addE_theta);
  addE_phi = new TH2F("Addback_Energy_vs_phi_segment", "Addback Energy vs. phi (segment)", 1800, -180, 180, 4096, 0, 8192);
  addE_phi->GetXaxis()->SetTitle("#phi (deg)");
  addE_phi->GetYaxis()->SetTitle("Addback Energy");
  tigList->Add(addE_phi);
  theta_phi = new TH2F("theta_vs_phi_segment", "theta vs. phi (segment, addback data)", 1800, 0, 180, 1800, -180, 180);
  theta_phi->GetXaxis()->SetTitle("#theta (deg)");
  theta_phi->GetYaxis()->SetTitle("#phi (deg)");
  tigList->Add(theta_phi);
  tigRate = new TH1F("TIGRESS_Total_Rate", "TIGRESS Total Rate", 8192, 0, 8192);
  tigRate->GetXaxis()->SetTitle("Run Time (s)");
  tigRate->GetYaxis()->SetTitle("Count/s");
  tigList->Add(tigRate);

  //Timing spectra
  tigT_tigT = new TH1F("Tigress_Tigress_time","Tigress - Tigress time",4096,-2048,2048); 
  tigT_tigT->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(tigT_tigT);
  addT_addT = new TH1F("Addback_Addback_time","Addback - Addback time",4096,-2048,2048); 
  addT_addT->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(addT_addT);
  tigT_tigT_tsep = new TH1F("Tigress_Tigress_time_tsep","Tigress - Tigress time, time separated",4096,-2048,2048); 
  tigT_tigT_tsep->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(tigT_tigT_tsep);
  addT_addT_tsep = new TH1F("Addback_Addback_time_tsep","Addback - Addback time, time separated",4096,-2048,2048); 
  addT_addT_tsep->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(addT_addT_tsep);
  tigT_tigT_tsepmult2 = new TH1F("Tigress_Tigress_time_tsepmult2","Tigress - Tigress time, time separated multiplicity 2",4096,-2048,2048); 
  tigT_tigT_tsepmult2->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(tigT_tigT_tsepmult2);
  addT_addT_tsepmult2 = new TH1F("Addback_Addback_time_tsepmult2","Addback - Addback time, time separated multiplicity 2",4096,-2048,2048); 
  addT_addT_tsepmult2->GetXaxis()->SetTitle("t_{TIGRESS} - t_{TIGRESS} (ns)");
  timingList->Add(addT_addT_tsepmult2);

  //TIGRESS-TIGRESS
  tigPos_tigPos = new TH2I("TIGRESS_TIGRESS_position", "TIGRESS-TIGRESS position", 64, 0, 64, 64, 0, 64);
  tigPos_tigPos->GetXaxis()->SetTitle("TIGRESS position 1");
  tigPos_tigPos->GetYaxis()->SetTitle("TIGRESS position 2");
  tigtigList->Add(tigPos_tigPos);
  tigE_tigE = new TH2F("TIGRESS_Gamma_Gamma", "TIGRESS Gamma-Gamma", 4096, 0, 4096, 4096, 0, 4096);
  tigE_tigE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  tigE_tigE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(tigE_tigE);
  addE_addE = new TH2F("Addback_Gamma_Gamma", "Addback Gamma-Gamma", 4096, 0, 4096, 4096, 0, 4096);
  addE_addE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(addE_addE);
  tigE_tigE_tsep = new TH2F("TIGRESS_Gamma_Gamma_tsep", "TIGRESS Gamma-Gamma time separated", 4096, 0, 4096, 4096, 0, 4096);
  tigE_tigE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  tigE_tigE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(tigE_tigE_tsep);
  addE_addE_tsep = new TH2F("Addback_Gamma_Gamma_tsep", "Addback Gamma-Gamma time separated", 4096, 0, 4096, 4096, 0, 4096);
  addE_addE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(addE_addE_tsep);
  tigE_tigE_tsepmult2 = new TH2F("TIGRESS_Gamma_Gamma_tsepmult2", "TIGRESS Gamma-Gamma time separated multiplicity 2", 4096, 0, 4096, 4096, 0, 4096);
  tigE_tigE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  tigE_tigE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(tigE_tigE_tsepmult2);
  addE_addE_tsepmult2 = new TH2F("Addback_Gamma_Gamma_tsepmult2", "Addback Gamma-Gamma time separated multiplicity 2", 4096, 0, 4096, 4096, 0, 4096);
  addE_addE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(addE_addE_tsepmult2);
  tigE_tigE_trand = new TH2F("TIGRESS_Gamma_Gamma_trand", "TIGRESS Gamma-Gamma time random", 4096, 0, 4096, 4096, 0, 4096);
  tigE_tigE_trand->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  tigE_tigE_trand->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(tigE_tigE_trand);
  addE_addE_trand = new TH2F("Addback_Gamma_Gamma_trand", "Addback Gamma-Gamma time random", 4096, 0, 4096, 4096, 0, 4096);
  addE_addE_trand->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  addE_addE_trand->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  tigtigList->Add(addE_addE_trand);

}
