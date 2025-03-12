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
TH2F *hpgeE_time, *addE_time;

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
    hpgeE_time = new TH2F("Gamma_Energy_vs_time", "Gamma Energy vs. time", 2048, 0, 2048, 32768, 0, 32768);
    hpgeE_time->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time);
    addE_time = new TH2F("Addback_Energy_vs_time", "Addback Energy vs. time", 2048, 0, 2048, 32768, 0, 32768);
    addE_time->GetXaxis()->SetTitle("Energy (keV)");
    addE_time->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(addE_time);


}
