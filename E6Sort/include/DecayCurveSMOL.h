#ifndef DecayCurveS_h
#define DecayCurveS_h

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

TList *dcList; 

//decay curve
TH1F *counts_time;
TH2F *hpgeE_time;

class DecayCurveS {

	public :

		DecayCurveS(){;} 
		void SortData(const char *, const uint64_t);
		void Initialise();
};
#endif

void DecayCurveS::Initialise() {

    cout << "Creating lists" << endl;

    dcList = new TList;

    cout << "Creating histograms" << endl;
    //histogram names shouldn't use spaces, to aid GRSISort-based analysis

    //decay curve Spectra
    counts_time = new TH1F("Counts_vs_time", "Gamma Energy vs. time", 32768, 0, 4096);
    counts_time->GetXaxis()->SetTitle("Time (min)");
    counts_time->GetYaxis()->SetTitle("Counts / min");
    dcList->Add(counts_time);
    hpgeE_time = new TH2F("Gamma_Energy_vs_time", "Gamma Energy vs. time", 2048, 0, 4096, 32768, 0, 4096);
    hpgeE_time->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time);


}
