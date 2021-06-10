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
#include "Declarations.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;

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
  tipPIDGatedList = new TList;
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
  num_addr = new TH2F("TIGRESS position vs. address", "TIGRESS position vs. address", 64, 0, 64, 16384, 0, 16384);
  tigList->Add(num_addr);
  addE_ANum_Downstream = new TH1F("Addback Energy for Downstream Lampshape positions","Addback Energy for Downstream Lampshape positions;Energy (keV);Counts/keV)",8192,0,8192);
  tigList->Add(addE_ANum_Downstream);
  addE_ANum_Upstream = new TH1F("Addback Energy for Upstream Lampshape positions","Addback Energy for Upstream Lampshape positions;Energy (keV);Counts/keV)",8192,0,8192);
  tigList->Add(addE_ANum_Upstream);
  addE_ANum_Corona = new TH1F("Addback Energy for Corona positions","Addback Energy for Corona positions;Energy (keV);Counts/keV)",8192,0,8192);
  tigList->Add(addE_ANum_Corona);
  tigTigDetectors = new TH2I("TIGRESS multiplicity 2 Array Number map","TIGRESS multiplicity 2 Array Number map;Hit 0 Array Number;Hit 1 Array Number",64, 0, 64,64, 0, 64);
  tigList->Add(tigTigDetectors);
  tigTigCloverDetectors = new TH2I("TIGRESS multiplicity 2 Clover Number map","TIGRESS multiplicity 2 Clover Number map;Hit 0 Clover Number;Hit 1 Clover Number",16, 0, 16,16, 0, 16);
  tigList->Add(tigTigCloverDetectors);

  //Raw TIP Spectra
  tip_pos = new TH1I("TIP position","TIP position",129,0,129); 
  tipList->Add(tip_pos);
  tip_ring = new TH1I("TIP ring","TIP ring",10,0,10);
  tipList->Add(tip_ring);
  tip_E = new TH1F("TIP energy","TIP energy",8192,0,8192); 
  tipList->Add(tip_E);
  tip_mult = new TH1I("TIP multiplicity","TIP multiplicity",10,0,10); 
  tipList->Add(tip_mult);
  tip_wfrmsize = new TH1F("TIP waveform size","TIP waveform size",4096,0,4096); 
  tipList->Add(tip_wfrmsize);
  tip_CFDFitDiff = new TH1F("TIP Fit - CFD timing","TIP Fit - CFD timing",2048,-1024,1024);
  tip_CFDFitDiff->GetXaxis()->SetTitle("t_{fit} - t_{CFD} (ns)");
  tipList->Add(tip_CFDFitDiff);
  tip_CFDFitDiffHighE = new TH1F("TIP Fit - CFD timing (E threshold above noise)","TIP Fit - CFD timing (E threshold above noise)",2048,-1024,1024);
  tip_CFDFitDiffHighE->GetXaxis()->SetTitle("t_{fit} - t_{CFD} (ns)");
  tipList->Add(tip_CFDFitDiffHighE);
  tip_E_pos = new TH2F("TIP energy vs position","TIP energy vs position",8192,0,8192,129,0,129);
  tip_E_pos->GetXaxis()->SetTitle("E (arb.)");
  tip_E_pos->GetYaxis()->SetTitle("TIP position");
  tipList->Add(tip_E_pos);
  tip_E_waveformAmp = new TH2F("TIP energy vs wavform amplitude","TIP energy vs wavform amplitude",4096,0,8192,4096,0,8192);
  tip_E_waveformAmp->GetXaxis()->SetTitle("E (arb.)");
  tip_E_waveformAmp->GetYaxis()->SetTitle("Waveform Amplitude (arb.)");
  tipList->Add(tip_E_waveformAmp);
  tiptipT = new TH1F("TIP-TIP timing","TIP-TIP timing",8192,-4096,4096);
  tipList->Add(tiptipT);
  tip_lastEvtTDiff_pos = new TH2F("TIP event time diff vs position","TIP event time diff vs position",8192,0,16384,129,0,129);
  tip_lastEvtTDiff_pos->GetXaxis()->SetTitle("Event time diff (us)");
  tip_lastEvtTDiff_pos->GetYaxis()->SetTitle("TIP position");
  tipList->Add(tip_lastEvtTDiff_pos);
  tipTipDetectors = new TH2I("TIP multiplicity 2 detector map","TIP multiplicity 2 detector map;Hit 0 Channel;Hit 1 Channel",129,0,129,129,0,129);
  tipList->Add(tipTipDetectors);
  tipPileupVsRing = new TH1F("TIP pileup fraction vs ring","TIP pileup fraction vs ring",10,0,10);
  tipPileupVsRing->GetXaxis()->SetTitle("TIP ring");
  tipPileupVsRing->GetYaxis()->SetTitle("Pileup fraction");
  tipList->Add(tipPileupVsRing);


  //TIP PID
  tip_E_PID_Sum = new TH2F("TIP energy vs PID (Sum)","TIP energy vs PID (Sum)",8192,0,8192,512,0,512);
  tip_E_PID_Sum->GetYaxis()->SetTitle("A_{S}/A_{F} x 100 + 100");
  tipPIDList->Add(tip_E_PID_Sum);
  for(int i=0; i<NTIP; i++){
    tip_E_PID[i] = new TH2F(Form("TIP energy vs PID (Pos %i)",i+1),Form("TIP energy vs PID (Pos %i)",i+1),8192,0,8192,512,0,512);
    tip_E_PID[i]->GetYaxis()->SetTitle("A_{S}/A_{F} x 100 + 100");
    tipPIDList->Add(tip_E_PID[i]);
  }

  //TIP PID Gated
  tip_particle_mult = new TH2I("TIP particle multiplicity","TIP particle multiplicity",10,0,10,10,0,10);
  tip_particle_mult->GetXaxis()->SetTitle("Proton multiplicity");
  tip_particle_mult->GetYaxis()->SetTitle("Alpha multiplicity");
  tipPIDGatedList->Add(tip_particle_mult);

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
  tigT_bgoT_supp= new TH1F("Tigress-BGO time (suppressed)","Tigress-BGO time (suppressed)",4096,-2048,2048);
  tigT_bgoT_supp->GetXaxis()->SetTitle("Same-detector Tigress-BGO time (ns)");
  tigbgoList->Add(tigT_bgoT_supp);

  //TIGRESS-TIP
  tiptig_mult = new TH2I("TIP-TIGRESS multiplicity","TIP-TIGRESS multiplicity",10,0,10,10,0,10);
  tiptig_mult->GetXaxis()->SetTitle("TIP multiplicity");
  tiptig_mult->GetYaxis()->SetTitle("TIGRESS multiplicity (suppressed)");
  tiptigList->Add(tiptig_mult);
  tipT_tigT_diff = new TH1F("TIP fit - Tigress time","TIP fit - Tigress time", 4096,-4096,4096);
  tipT_tigT_diff->GetXaxis()->SetTitle("t_{TIP, fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipT_tigT_diff);
  tipTCFD_tigT_diff = new TH1F("TIP CFD - Tigress time","TIP CFD - Tigress time", 4096,-4096,4096);
  tipTCFD_tigT_diff->GetXaxis()->SetTitle("t_{TIP, CFD} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipTCFD_tigT_diff);

  tigE_TIPtg = new TH1F("Tigress energy (CFD time gated)","Tigress energy (CFD time gated)", 8192, 0, 8192); 
  tiptigList->Add(tigE_TIPtg);
  tipE_tigtg = new TH1F("TIP energy (CFD time gated)","TIP energy (CFD time gated)", 8192, 0, 8192); 
  tiptigList->Add(tipE_tigtg);
  tipT_tigT_difftg = new TH1F("TIP fit - Tigress time (CFD time gated)","TIP fit - Tigress time (CFD time gated)", 4096,-4096,4096);
  tipT_tigT_difftg->GetXaxis()->SetTitle("t_{TIP,fit} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipT_tigT_difftg);
  tipTCFD_tigT_difftg = new TH1F("TIP CFD - Tigress time (CFD time gated)","TIP fit - Tigress time (CFD time gated)", 4096,-4096,4096);
  tipTCFD_tigT_difftg->GetXaxis()->SetTitle("t_{TIP,CFD} - t_{TIGRESS} (ns)");
  tiptigList->Add(tipTCFD_tigT_difftg);

  tigE_0p2aGate = new TH1F("Tigress energy (0p2a gated)","Tigress energy (0p2a gated)", 8192, 0, 8192);
  tiptigList->Add(tigE_0p2aGate);
  tigE_2p0aGate = new TH1F("Tigress energy (2p0a gated)","Tigress energy (2p0a gated)", 8192, 0, 8192);
  tiptigList->Add(tigE_2p0aGate);

  tigE_0p2aDetSumGate = new TH1F("Tigress energy (0p2a - rough - common PID gate on all TIP positions)","Tigress energy (0p2a - rough - common PID gate on all TIP positions)", 8192, 0, 8192);
  tiptigList->Add(tigE_0p2aDetSumGate);
  tigE_2p0aDetSumGate = new TH1F("Tigress energy (2p0a - rough - common PID gate on all TIP positions)","Tigress energy (2p0a - rough - common PID gate on all TIP positions)", 8192, 0, 8192);
  tiptigList->Add(tigE_2p0aDetSumGate);
  tigE_horDetSumGate = new TH1F("Tigress energy (horizontal gate - rough - common PID Gate on all TIP positions)","Tigress energy (horizontal gate - rough - common PID Gate on all TIP positions);#gamma Energy (keV);Counts/keV",8192,0,8192);
  tiptigList->Add(tigE_horDetSumGate);
  tipChannel_HorizontalCut = new TH1F("TIP Channel - Horizontal Cut","TIP Channel - Horizontal Cut;Channel Number;",129,0,129);
  tiptigList->Add(tipChannel_HorizontalCut);

  tigE_2p0aDetSumGate_Downstream = new TH1F("Downstream Tigress energy (2p0a - rough - common PID gate on all TIP positions)","Downstream Tigress energy (2p0a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_2p0aDetSumGate_Downstream);
  tigE_2p0aDetSumGate_Upstream = new TH1F("Upstream Tigress energy (2p0a - rough - common PID gate on all TIP positions)","Upstream Tigress energy (2p0a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_2p0aDetSumGate_Upstream);
  tigE_2p0aDetSumGate_Corona = new TH1F("Corona Tigress energy (2p0a - rough - common PID gate on all TIP positions)","Corona Tigress energy (2p0a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_2p0aDetSumGate_Corona);

  tigE_0p2aDetSumGate_Downstream = new TH1F("Downstream Tigress energy (0p2a - rough - common PID gate on all TIP positions)","Downstream Tigress energy (0p2a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_0p2aDetSumGate_Downstream);
  tigE_0p2aDetSumGate_Upstream = new TH1F("Upstream Tigress energy (0p2a - rough - common PID gate on all TIP positions)","Upstream Tigress energy (0p2a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_0p2aDetSumGate_Upstream);
  tigE_0p2aDetSumGate_Corona = new TH1F("Corona Tigress energy (0p2a - rough - common PID gate on all TIP positions)","Corona Tigress energy (0p2a - rough - common PID gate on all TIP positions);Energy (keV);Counts/keV", 8192, 0, 8192);
  tiptigList->Add(tigE_0p2aDetSumGate_Corona);


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

  //common gate for all detectors
  /* //old cuts that are backwards (alpha cut is actually on proton channel)
  alphaPosCutSum = new TCutG("common alpha cut",6);
  alphaPosCutSum->SetPoint(0,225,147);
  alphaPosCutSum->SetPoint(1,274,241);
  alphaPosCutSum->SetPoint(2,3675,320);
  alphaPosCutSum->SetPoint(3,4027,203);
  alphaPosCutSum->SetPoint(4,596,144);
  alphaPosCutSum->SetPoint(5,225,147);
  protonPosCutSum = new TCutG("common proton cut",5);
  protonPosCutSum->SetPoint(0,1022,111);
  protonPosCutSum->SetPoint(1,905,145);
  protonPosCutSum->SetPoint(2,5258,199);
  protonPosCutSum->SetPoint(3,4720,118);
  protonPosCutSum->SetPoint(4,1022,111);
  */

  /*//cuts from run sum 53951-53966

  alphaPosCutSum = new TCutG("common alpha cut",6);
  alphaPosCutSum->SetPoint(0,267,96);
  alphaPosCutSum->SetPoint(1,3986,143);
  alphaPosCutSum->SetPoint(2,4078,164);
  alphaPosCutSum->SetPoint(3,120,137);
  alphaPosCutSum->SetPoint(4,167,99);
  alphaPosCutSum->SetPoint(5,267,96);

  protonPosCutSum = new TCutG("common proton cut",6);
  protonPosCutSum->SetPoint(0,578,149);
  protonPosCutSum->SetPoint(1,1952,216);
  protonPosCutSum->SetPoint(2,1806,250);
  protonPosCutSum->SetPoint(3,305,202);
  protonPosCutSum->SetPoint(4,259,166);
  protonPosCutSum->SetPoint(5,578,149);*/

  //cuts from run sum 53996-54004

  alphaPosCutSum = new TCutG("common alpha cut",6);
  alphaPosCutSum->SetPoint(0,267,96);
  alphaPosCutSum->SetPoint(1,3986,143);
  alphaPosCutSum->SetPoint(2,4078,164);
  alphaPosCutSum->SetPoint(3,120,137);
  alphaPosCutSum->SetPoint(4,167,99);
  alphaPosCutSum->SetPoint(5,267,96);

  protonPosCutSum = new TCutG("common proton cut",6);
  protonPosCutSum->SetPoint(0,596,141);
  protonPosCutSum->SetPoint(1,2673,179);
  protonPosCutSum->SetPoint(2,2383,258);
  protonPosCutSum->SetPoint(3,305,202);
  protonPosCutSum->SetPoint(4,237,139);
  protonPosCutSum->SetPoint(5,596,141);


  tipPIDGateList->Add(alphaPosCutSum);
  tipPIDGateList->Add(protonPosCutSum);


  //cut on strange horizontal feature
  horPosCutSum = new TCutG("common horizontal cut",8);
  horPosCutSum->SetPoint(0,213,45);
  horPosCutSum->SetPoint(1,8065,45);
  horPosCutSum->SetPoint(2,8060,58);
  horPosCutSum->SetPoint(3,843,67);
  horPosCutSum->SetPoint(4,466,81);
  horPosCutSum->SetPoint(5,194,72);
  horPosCutSum->SetPoint(6,188,64);
  horPosCutSum->SetPoint(7,213,45);

  tipPIDGateList->Add(horPosCutSum);


  //individual gates for each detector
  for(int i=0; i<NTIP; i++){
    alphaPosCut[i] = new TCutG(Form("pos %i alpha cut",i+1),7);
    protonPosCut[i] = new TCutG(Form("pos %i proton cut",i+1),7);

    //CHANGE LATER!! this is just proof of concept
    alphaPosCut[i]->SetPoint(0,4360,52);
    alphaPosCut[i]->SetPoint(1,4399,150);
    alphaPosCut[i]->SetPoint(2,6513,197);
    alphaPosCut[i]->SetPoint(3,6867,51);
    alphaPosCut[i]->SetPoint(4,5901,56);
    alphaPosCut[i]->SetPoint(5,5066,27);
    alphaPosCut[i]->SetPoint(6,4360,52);

    protonPosCut[i]->SetPoint(0,4360,52);
    protonPosCut[i]->SetPoint(1,4399,150);
    protonPosCut[i]->SetPoint(2,6513,197);
    protonPosCut[i]->SetPoint(3,6867,51);
    protonPosCut[i]->SetPoint(4,5901,56);
    protonPosCut[i]->SetPoint(5,5066,27);
    protonPosCut[i]->SetPoint(6,4360,52);
  }

  //manual gate definitions go here

  //add gates to list
  for(int i=0; i<NTIP; i++){
    tipPIDGateList->Add(alphaPosCut[i]);
    tipPIDGateList->Add(protonPosCut[i]);
  }


}
