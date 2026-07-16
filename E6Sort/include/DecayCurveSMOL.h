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

double addbackE[NGRIFPOS],maxABHitE[NGRIFPOS];
double addbackT[NGRIFPOS];

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree separation
uint8_t hitMap180degAB[NGRIFPOS][NGRIFPOS]; //1st index = clover of hit, 2nd index = clover of 2nd hit, val = 1 indicates 180 degree separation

TApplication *theApp;

TList *dcList, *dcABList; 

//decay curve
TH1D *counts_time, *counts_time_1477, *counts_time_934, *counts_time_934_bg, *counts_time_685;
TH1D *counts_time_511betaplus, *counts_time_511betaplus_934, *counts_time_511betaplus_1238, *counts_time_511betaplusAB, *counts_time_511betaplusAB_934, *counts_time_511betaplusAB_1238, *counts_time_511betaplusAB_955;
TH2D *hpgeE_time, *addE_time, *angle511betaplus_time, *angle511betaplus_ecoinc, *angle511betaplusAB_time, *angle511betaplusAB_ecoinc;
TH2D *hpgeE_time_betaplus, *hpgeE_time_betaplus_934, *hpgeE_time_betaplus_3rdgamma, *addE_time_betaplus, *addE_time_betaplus_934, *addE_time_betaplus_955, *addE_time_betaplus_3rdgamma;
TH1D *hpgeE, *addE, *hpgeE934, *addE934;

#endif

void InitialiseHists() {

    cout << "Creating lists" << endl;

    dcList = new TList;
    dcABList = new TList;

    cout << "Creating histograms" << endl;
    //histogram names shouldn't use spaces, to aid GRSISort-based analysis

    //decay curve addback spectra
    addE = new TH1D("Addback_Energy", "Addback Energy", 16384, 0, 8192);
    addE->GetXaxis()->SetTitle("Energy (keV)");
    addE->GetYaxis()->SetTitle("Counts");
    dcABList->Add(addE);
    addE934 = new TH1D("Addback_Energy_934", "Addback Energy (934 keV gated)", 16384, 0, 8192);
    addE934->GetXaxis()->SetTitle("Energy (keV)");
    addE934->GetYaxis()->SetTitle("Counts");
    dcABList->Add(addE934);
    counts_time = new TH1D("Counts_vs_time", "Addback Energy vs. time", 32768, 0, 32768);
    counts_time->GetXaxis()->SetTitle("Time (min)");
    counts_time->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time);
    counts_time_685 = new TH1D("Counts_vs_time_685", "Addback Energy vs. time (685 keV gate)", 32768, 0, 32768);
    counts_time_685->GetXaxis()->SetTitle("Time (min)");
    counts_time_685->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_685);
    counts_time_934 = new TH1D("Counts_vs_time_934", "Addback Energy vs. time (934 keV gate)", 32768, 0, 32768);
    counts_time_934->GetXaxis()->SetTitle("Time (min)");
    counts_time_934->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_934);
    counts_time_934_bg = new TH1D("Counts_vs_time_934_bg", "Addback Energy vs. time (934 keV Compton background gate)", 32768, 0, 32768);
    counts_time_934_bg->GetXaxis()->SetTitle("Time (min)");
    counts_time_934_bg->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_934_bg);
    counts_time_1477 = new TH1D("Counts_vs_time_1477", "Addback Energy vs. time (1477 keV gate)", 32768, 0, 32768);
    counts_time_1477->GetXaxis()->SetTitle("Time (min)");
    counts_time_1477->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_1477);
    counts_time_511betaplusAB = new TH1D("Counts_vs_time_511betaplusAB", "Addback Energy vs. time (511 keV 180 degree coinc. from beta+, clover addback)", 32768, 0, 32768);
    counts_time_511betaplusAB->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplusAB->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_511betaplusAB);
    counts_time_511betaplusAB_1238 = new TH1D("Counts_vs_time_511betaplusAB_1238", "Addback Energy vs. time (511 keV 180 degree coinc. from beta+, clover addback, 1238 keV gate)", 32768, 0, 32768);
    counts_time_511betaplusAB_1238->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplusAB_1238->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_511betaplusAB_1238);
    counts_time_511betaplusAB_934 = new TH1D("Counts_vs_time_511betaplusAB_934", "Addback Energy vs. time (511 keV 180 degree coinc. from beta+, clover addback, 934 keV gate)", 32768, 0, 32768);
    counts_time_511betaplusAB_934->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplusAB_934->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_511betaplusAB_934);
    counts_time_511betaplusAB_955 = new TH1D("Counts_vs_time_511betaplusAB_955", "Addback Energy vs. time (511 keV 180 degree coinc. from beta+, clover addback, 955 keV gate)", 32768, 0, 32768);
    counts_time_511betaplusAB_955->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplusAB_955->GetYaxis()->SetTitle("Counts / minute");
    dcABList->Add(counts_time_511betaplusAB_955);
    addE_time = new TH2D("Addback_Energy_vs_time", "Addback Energy vs. time", 2048, 0, 4096, 32768, 0, 32768);
    addE_time->GetXaxis()->SetTitle("Energy (keV)");
    addE_time->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(addE_time);
    addE_time_betaplus = new TH2D("Addback_Energy_vs_time_betaplus", "Addback Energy vs. time (180 degree coinc. with 511 keV)", 2048, 0, 4096, 32768, 0, 32768);
    addE_time_betaplus->GetXaxis()->SetTitle("Energy (keV)");
    addE_time_betaplus->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(addE_time_betaplus);
    addE_time_betaplus_934 = new TH2D("Addback_Energy_vs_time_betaplus_934", "Addback Energy vs. time (180 degree coinc. with 511 keV, 934 keV gate)", 2048, 0, 4096, 32768, 0, 32768);
    addE_time_betaplus_934->GetXaxis()->SetTitle("Energy (keV)");
    addE_time_betaplus_934->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(addE_time_betaplus_934);
    addE_time_betaplus_955 = new TH2D("Addback_Energy_vs_time_betaplus_955", "Addback Energy vs. time (180 degree coinc. with 511 keV, 955 keV gate)", 2048, 0, 4096, 32768, 0, 32768);
    addE_time_betaplus_955->GetXaxis()->SetTitle("Energy (keV)");
    addE_time_betaplus_955->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(addE_time_betaplus_955);
    addE_time_betaplus_3rdgamma = new TH2D("Addback_Energy_vs_time_betaplus_3rdgamma", "Addback Energy vs. time (180 degree coinc. with 511 keV, 3rd gamma energy)", 2048, 0, 4096, 32768, 0, 32768);
    addE_time_betaplus_3rdgamma->GetXaxis()->SetTitle("Energy (keV)");
    addE_time_betaplus_3rdgamma->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(addE_time_betaplus_3rdgamma);
    angle511betaplusAB_ecoinc = new TH2D("angle511betaplusAB_ecoinc", "Addback 511 keV angle vs. coincident energy", 256, -256, 256, 16384, 0, 8192);
    angle511betaplusAB_ecoinc->GetXaxis()->SetTitle("Angle between 511 keV gammas (deg)");
    angle511betaplusAB_ecoinc->GetYaxis()->SetTitle("Coincident energy (keV)");
    dcABList->Add(angle511betaplusAB_ecoinc);
    angle511betaplusAB_time = new TH2D("angle511betaplusAB_time", "Addback 511 keV angle vs. time", 256, -256, 256, 32768, 0, 32768);
    angle511betaplusAB_time->GetXaxis()->SetTitle("Angle between 511 keV gammas (deg)");
    angle511betaplusAB_time->GetYaxis()->SetTitle("Time (min)");
    dcABList->Add(angle511betaplusAB_time);

    //decay curve non-addback (single crystal) spectra
    hpgeE = new TH1D("NonAddback_Energy", "Single-crystal Energy", 16384, 0, 8192);
    hpgeE->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE->GetYaxis()->SetTitle("Counts");
    dcList->Add(hpgeE);
    hpgeE934 = new TH1D("NonAddback_Energy_934", "Single-crystal Energy (934 keV gated)", 16384, 0, 8192);
    hpgeE934->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE934->GetYaxis()->SetTitle("Counts");
    dcList->Add(hpgeE934);
    counts_time_511betaplus = new TH1D("Counts_vs_time_511betaplus", "Addback Energy vs. time (511 keV 180 degree coinc. from beta, non-addback)", 32768, 0, 32768);
    counts_time_511betaplus->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplus->GetYaxis()->SetTitle("Counts / minute");
    dcList->Add(counts_time_511betaplus);
    counts_time_511betaplus_1238 = new TH1D("Counts_vs_time_511betaplus_1238", "Addback Energy vs. time (511 keV 180 degree coinc. from beta, non-addback, 1238 keV gate)", 32768, 0, 32768);
    counts_time_511betaplus_1238->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplus_1238->GetYaxis()->SetTitle("Counts / minute");
    dcList->Add(counts_time_511betaplus_1238);
    counts_time_511betaplus_934 = new TH1D("Counts_vs_time_511betaplus_934", "Addback Energy vs. time (511 keV 180 degree coinc. from beta, non-addback, 934 keV gate)", 32768, 0, 32768);
    counts_time_511betaplus_934->GetXaxis()->SetTitle("Time (min)");
    counts_time_511betaplus_934->GetYaxis()->SetTitle("Counts / minute");
    dcList->Add(counts_time_511betaplus_934);
    hpgeE_time = new TH2D("NonAddback_Energy_vs_time", "Single-crystal Energy vs. time", 2048, 0, 4096, 32768, 0, 32768);
    hpgeE_time->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time);
    hpgeE_time_betaplus = new TH2D("NonAddback_Energy_vs_time_betaplus", "Single-crystal Energy vs. time (180 degree coinc. with 511 keV)", 2048, 0, 4096, 32768, 0, 32768);
    hpgeE_time_betaplus->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time_betaplus->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time_betaplus);
    hpgeE_time_betaplus_934 = new TH2D("NonAddback_Energy_vs_time_betaplus_934", "Single-crystal Energy vs. time (180 degree coinc. with 511 keV, 934 keV gate)", 2048, 0, 4096, 32768, 0, 32768);
    hpgeE_time_betaplus_934->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time_betaplus_934->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time_betaplus_934);
    hpgeE_time_betaplus_3rdgamma = new TH2D("NonAddback_Energy_vs_time_betaplus_3rdgamma", "Single-crystal Energy vs. time (180 degree coinc. with 511 keV, 3rd gamma energy)", 2048, 0, 4096, 32768, 0, 32768);
    hpgeE_time_betaplus_3rdgamma->GetXaxis()->SetTitle("Energy (keV)");
    hpgeE_time_betaplus_3rdgamma->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(hpgeE_time_betaplus_3rdgamma);
    angle511betaplus_ecoinc = new TH2D("angle511betaplus_ecoinc", "511 keV angle vs. coincident energy", 256, -256, 256, 16384, 0, 8192);
    angle511betaplus_ecoinc->GetXaxis()->SetTitle("Angle between 511 keV gammas (deg)");
    angle511betaplus_ecoinc->GetYaxis()->SetTitle("Coincident energy (keV)");
    dcList->Add(angle511betaplus_ecoinc);
    angle511betaplus_time = new TH2D("angle511betaplus_time", "511 keV angle vs. time", 256, -256, 256, 32768, 0, 32768);
    angle511betaplus_time->GetXaxis()->SetTitle("Angle between 511 keV gammas (deg)");
    angle511betaplus_time->GetYaxis()->SetTitle("Time (min)");
    dcList->Add(angle511betaplus_time);

}
