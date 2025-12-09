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
TH1D *hpgeE;
TH2D *hpgeE_ANum;
TH1I *hpgeMult;

//Timing
TH1I *hpgeT_hpgeT_le;
TH1D *hpgeT_hpgeT;
TH1D *hpgeT_hpgeT_tsep;
TH1D *hpgeT_hpgeT_tsepmult2;
TH1D *hpgeT_hpgeT_tseprand;
TH2D *hpgeT_hpgeT_EDiff, *hpgeT_hpgeT_hpgeE, *hpgeT_hpgeT_hpgeE_NoCFDfail, *hpgeT_hpgeT_hpgeE_1CFDfail, *hpgeT_hpgeT_hpgeE_2CFDfail;

//HPGe-HPGe
TH1D *hpge_hpge_dist, *hpge_hpge_angle;
TH2I *hpgePos_hpgePos, *hpgePos_hpgePos_lowEDiff;
TH2D *hpgeE_hpgeE;
TH2D *hpgeE_hpgeE_tsep;
TH2D *hpgeE_hpgeE_tsepmult2;
TH2D *hpgeE_hpgeE_tseprand;
TH2D *hpgeE_hpgeE_180deg, *hpgeE_hpgeE_180deg_1477gate;
TH1D *hpgeE_hpgeE_180deg_proj;
TH1D *hpgeE_hpgeE_180deg_sum, *hpgeE_hpgeE_180deg_sum_1477gate;
TH2D *hpgeE_hpgeE_180deg_sum_tDiff, *hpgeE_hpgeE_180deg_sum_tDiff_1477gate, *hpgeE_hpgeE_180deg_sum_tDiff_685gate;

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
  hpgeE = new TH1D("HPGe_Energy", "HPGe Energy (non-addback)", 16384, 0, 8192);
  hpgeList->Add(hpgeE);
  hpgeE_ANum = new TH2D("HPGe_Energy_vs_Array_Number", "HPGe Energy vs. Array Number", 64, 0, 64, 8192, 0, 8192);
  hpgeList->Add(hpgeE_ANum);
  hpgeMult = new TH1I("hpgeMult", "HPGe multiplicity (non-addback)", 64, 0, 64);
  hpgeMult->GetYaxis()->SetTitle("Counts");
  hpgeList->Add(hpgeMult);

  //Timing spectra
  hpgeT_hpgeT_le = new TH1I("HPGe_HPGe_tstime","HPGe - HPGe leading edge (timestamp) time",512,-256,256); 
  hpgeT_hpgeT_le->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (timestamp units)");
  timingList->Add(hpgeT_hpgeT_le);
  hpgeT_hpgeT = new TH1D("HPGe_HPGe_time","HPGe - HPGe time",4096,-2048,2048); 
  hpgeT_hpgeT->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT);
  hpgeT_hpgeT_tsep = new TH1D("HPGe_HPGe_time_tsep","HPGe - HPGe time, time separated",4096,-2048,2048); 
  hpgeT_hpgeT_tsep->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsep);
  hpgeT_hpgeT_tsepmult2 = new TH1D("HPGe_HPGe_time_tsepmult2","HPGe - HPGe time, time separated multiplicity 2",4096,-2048,2048); 
  hpgeT_hpgeT_tsepmult2->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tsepmult2);
  hpgeT_hpgeT_tseprand = new TH1D("HPGe_HPGe_time_tseprand","HPGe - HPGe time, time-random separated",4096,-2048,2048); 
  hpgeT_hpgeT_tseprand->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  timingList->Add(hpgeT_hpgeT_tseprand);
  hpgeT_hpgeT_EDiff = new TH2D("HPGe_HPGe_time_vs_e_diff", "HPGe - HPGe time vs. energy difference", 4096, -2048, 2048, 4096, -256, 256);
  hpgeT_hpgeT_EDiff->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  hpgeT_hpgeT_EDiff->GetYaxis()->SetTitle("E_{#gamma 2} - E_{#gamma 1} (keV)");
  timingList->Add(hpgeT_hpgeT_EDiff);
  hpgeT_hpgeT_hpgeE = new TH2D("HPGe_HPGe_time_vs_energy", "HPGe - HPGe time vs. energy", 4096, -2048, 2048, 4096, 0, 8192);
  hpgeT_hpgeT_hpgeE->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  hpgeT_hpgeT_hpgeE->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  timingList->Add(hpgeT_hpgeT_hpgeE);
  hpgeT_hpgeT_hpgeE_NoCFDfail = new TH2D("HPGe_HPGe_time_vs_energy_NoCFDfail", "HPGe - HPGe time vs. energy (no CFD fail)", 4096, -2048, 2048, 2048, 0, 4096);
  hpgeT_hpgeT_hpgeE_NoCFDfail->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  hpgeT_hpgeT_hpgeE_NoCFDfail->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  timingList->Add(hpgeT_hpgeT_hpgeE_NoCFDfail);
  hpgeT_hpgeT_hpgeE_1CFDfail = new TH2D("HPGe_HPGe_time_vs_energy_1CFDfail", "HPGe - HPGe time vs. energy (1 CFD fail hit)", 4096, -2048, 2048, 2048, 0, 4096);
  hpgeT_hpgeT_hpgeE_1CFDfail->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  hpgeT_hpgeT_hpgeE_1CFDfail->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  timingList->Add(hpgeT_hpgeT_hpgeE_1CFDfail);
  hpgeT_hpgeT_hpgeE_2CFDfail = new TH2D("HPGe_HPGe_time_vs_energy_2CFDfail", "HPGe - HPGe time vs. energy (2 CFD fail hits)", 4096, -2048, 2048, 2048, 0, 4096);
  hpgeT_hpgeT_hpgeE_2CFDfail->GetXaxis()->SetTitle("t_{HPGe} - t_{HPGe} (ns)");
  hpgeT_hpgeT_hpgeE_2CFDfail->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  timingList->Add(hpgeT_hpgeT_hpgeE_2CFDfail);

  //HPGe-HPGe
  hpge_hpge_dist = new TH1D("HPGe_HPGe_distance","HPGe-HPGe distance",512,0,512);
  hpge_hpge_dist->GetXaxis()->SetTitle("HPGe - HPGe distance (mm)");
  hpgehpgeList->Add(hpge_hpge_dist);
  hpge_hpge_angle = new TH1D("HPGe_HPGe_angle","HPGe-HPGe angle",512,0,360);
  hpge_hpge_angle->GetXaxis()->SetTitle("HPGe - HPGe angle (deg)");
  hpgehpgeList->Add(hpge_hpge_angle);
  hpgePos_hpgePos = new TH2I("HPGe_HPGe_position", "HPGe-HPGe position", 64, 0, 64, 64, 0, 64);
  hpgePos_hpgePos->GetXaxis()->SetTitle("HPGe position 1");
  hpgePos_hpgePos->GetYaxis()->SetTitle("HPGe position 2");
  hpgehpgeList->Add(hpgePos_hpgePos);
  hpgePos_hpgePos_lowEDiff = new TH2I("HPGe_HPGe_position_lowEDiff", "HPGe-HPGe position, Ediff < 0.1 keV", 64, 0, 64, 64, 0, 64);
  hpgePos_hpgePos_lowEDiff->GetXaxis()->SetTitle("HPGe position 1");
  hpgePos_hpgePos_lowEDiff->GetYaxis()->SetTitle("HPGe position 2");
  hpgehpgeList->Add(hpgePos_hpgePos_lowEDiff);
  hpgeE_hpgeE = new TH2D("HPGe_Gamma_Gamma", "HPGe Gamma-Gamma", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE);
  hpgeE_hpgeE_tsep = new TH2D("HPGe_Gamma_Gamma_tsep", "HPGe Gamma-Gamma time separated", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tsep->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsep->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsep);
  hpgeE_hpgeE_tsepmult2 = new TH2D("HPGe_Gamma_Gamma_tsepmult2", "HPGe Gamma-Gamma time separated multiplicity 2", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tsepmult2->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tsepmult2->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tsepmult2);
  hpgeE_hpgeE_tseprand = new TH2D("HPGe_Gamma_Gamma_tseprand", "HPGe Gamma-Gamma time-random separated", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_tseprand->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_tseprand->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_tseprand);
  hpgeE_hpgeE_180deg = new TH2D("HPGe_Gamma_Gamma_180deg", "HPGe Gamma-Gamma 180 degree coincidences", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_180deg->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_180deg->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg);
  hpgeE_hpgeE_180deg_1477gate = new TH2D("HPGe_Gamma_Gamma_180deg_1477gate", "HPGe Gamma-Gamma 180 degree coincidences, 1477 keV gate", 8192, 0, 4096, 8192, 0, 4096);
  hpgeE_hpgeE_180deg_1477gate->GetXaxis()->SetTitle("E_{#gamma 1} (keV)");
  hpgeE_hpgeE_180deg_1477gate->GetYaxis()->SetTitle("E_{#gamma 2} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_1477gate);
  hpgeE_hpgeE_180deg_proj = new TH1D("HPGe_Gamma_Gamma_180deg_proj", "HPGe Gamma-Gamma 180 degree coincidences (projection)", 8192, 0, 4096);
  hpgeE_hpgeE_180deg_proj->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_proj);
  hpgeE_hpgeE_180deg_sum = new TH1D("HPGe_Gamma_Gamma_180deg_sum", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits)", 8192, 0, 4096);
  hpgeE_hpgeE_180deg_sum->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum);
  hpgeE_hpgeE_180deg_sum_1477gate = new TH1D("HPGe_Gamma_Gamma_180deg_sum_1477gate", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits), 1477 keV gate", 8192, 0, 4096);
  hpgeE_hpgeE_180deg_sum_1477gate->GetXaxis()->SetTitle("E_{#gamma} (keV)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum_1477gate);
  hpgeE_hpgeE_180deg_sum_tDiff = new TH2D("HPGe_Gamma_Gamma_180deg_sum_tDiff", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits vs. tDiff)", 8192, 0, 4096, 512,-256,256);
  hpgeE_hpgeE_180deg_sum_tDiff->GetXaxis()->SetTitle("E_{#gamma } (keV)");
  hpgeE_hpgeE_180deg_sum_tDiff->GetYaxis()->SetTitle("t_{HPGe} - t_{HPGe} (timestamp units)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum_tDiff);
  hpgeE_hpgeE_180deg_sum_tDiff_1477gate = new TH2D("HPGe_Gamma_Gamma_180deg_sum_tDiff_1477gate", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits vs. tDiff), 1477 keV gate", 8192, 0, 4096, 512,-256,256);
  hpgeE_hpgeE_180deg_sum_tDiff_1477gate->GetXaxis()->SetTitle("E_{#gamma } (keV)");
  hpgeE_hpgeE_180deg_sum_tDiff_1477gate->GetYaxis()->SetTitle("t_{HPGe} - t_{HPGe} (timestamp units)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum_tDiff_1477gate);
  hpgeE_hpgeE_180deg_sum_tDiff_685gate = new TH2D("HPGe_Gamma_Gamma_180deg_sum_tDiff_685gate", "HPGe Gamma-Gamma 180 degree coincidences (sum of hits vs. tDiff), 685 keV gate", 8192, 0, 4096, 512,-256,256);
  hpgeE_hpgeE_180deg_sum_tDiff_685gate->GetXaxis()->SetTitle("E_{#gamma } (keV)");
  hpgeE_hpgeE_180deg_sum_tDiff_685gate->GetYaxis()->SetTitle("t_{HPGe} - t_{HPGe} (timestamp units)");
  hpgehpgeList->Add(hpgeE_hpgeE_180deg_sum_tDiff_685gate);

}
