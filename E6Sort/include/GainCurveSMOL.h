#ifndef GainCurveS_h
#define GainCurveS_h

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

TList *gcList; 

//gain curve
TH2F *hpgeE_time[NTIGPOS*4];

class GainCurveS {

	public :

		GainCurveS(){;} 
		void SortData(const char *, const uint64_t);
		void Initialise();
};
#endif

void GainCurveS::Initialise() {

    cout << "Creating lists" << endl;

    gcList = new TList;

    cout << "Creating histograms" << endl;
    //histogram names shouldn't use spaces, to aid GRSISort-based analysis

    //gain curve Spectra
    for(int i=0; i<(NTIGPOS*4); i++){
        hpgeE_time[i] = new TH2F(Form("Gamma_Energy_time_%i",i), Form("Gamma Energy vs. time (crystal %i)",i), 8192, 0, 32768, 1024, 0, 1024);
        hpgeE_time[i]->GetXaxis()->SetTitle("Time (min)");
        hpgeE_time[i]->GetYaxis()->SetTitle("Energy (keV)");
        gcList->Add(hpgeE_time[i]);
    }


}
