#ifndef SortCode_h
#define SortCode_h

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TEmma.h"
#include "TParserLibrary.h"
#include "TEnv.h"
using namespace std;

TList *multList, *emmaList, *icList, *tigList, *addList, *tigtigList, *tigbgoList, *tigemmaList, *addemmaList, *timegatedList, *massgatedList, *ggatedList, *ggatedicList, *ggatedICSiList; 

TH2F *pgac, *icN, *iC0V1, *iC0V2, *iC0V3, *iC1V2, *iC1V3, *iC2V3, *tDE, *xPosT, *yPosT, *siET;
TH1F *xPos, *yPos, *siE, *icSum; 
TH2F *icSumVSi, *ic0VSi, *ic1VSi, *ic2VSi, *ic3VSi;
TH1F *icE[5], *ssbE[2];
TH2F *ssbET[2];

TH1F *tigE, *tigDop;
TH2F *tigE_ANum, *tigDop_ANum, *tigE_theta, *tigDop_theta; 
//begin CRN
TH2F* tigE_tigDop;
//end CRN

TH1F *addE, *addDop;
TH2F *addE_ANum, *addDop_ANum, *addE_theta, *addDop_theta, *tigDop_ANum_if_emma; 

TH1F *tigemmatof, *addemmatof;
TH2F *tigE_tof, *addE_tof, *tigDop_tof, *addDop_tof;

TH1F *tigE_tg,*addDop_tg,*addDop_tg_bg,*addDop_tg_bgs;
TH2F *addDopp_addDopp_tg, *tigE_xPos_tg, *tigDop_xPos_tg, *tigE_xPos_bg, *tigDop_xPos_bg, *addE_xPos_tg, *addDop_xPos_tg, *addE_xPos_bg, *addDop_xPos_bg, *pgac_tg, *pgac_bg, *icN_tg, *icN_bg;

TH1F *tigE_massG[2], *tigDop_massG[2], *addE_massG[2], *addDop_massG[2], *tigE_massG_bg[2], *tigDop_massG_bg[2], *addE_massG_bg[2], *addDop_massG_bg[2], *addDop_mtof[2];

TH2F *pgac_gGate[5], *pgac_aGate[5], *pgac_gGate_bg[5], *pgac_aGate_bg[5], *icN_aGate_tg[5], *icN_aGate_bg[5], *iC0V1_aGate_tg[5], *iC0V2_aGate_tg[5], *iC0V3_aGate_tg[5], *iC1V2_aGate_tg[5], *iC1V3_aGate_tg[5], *iC2V3_aGate_tg[5], *iC0V1_aGate_bg[5], *iC0V2_aGate_bg[5], *iC0V3_aGate_bg[5], *iC1V2_aGate_bg[5], *iC1V3_aGate_bg[5], *iC2V3_aGate_bg[5];

TH2F * icSumVSi_aGate_tg[5], *ic0VSi_aGate_tg[5], *ic1VSi_aGate_tg[5], *ic2VSi_aGate_tg[5], *ic3VSi_aGate_tg[5], * icSumVSi_aGate_bg[5], *ic0VSi_aGate_bg[5], *ic1VSi_aGate_bg[5], *ic2VSi_aGate_bg[5], *ic3VSi_aGate_bg[5];
TCutG *massGate[2], *icGate[2];

//TIGRESS-TIGRESS
TH1F *addT_addT; 
TH2F *addE_addE;

//TIGRESS-BGO
TH1I *bgo_mult, *bgo_det;
TH1F *tigT_bgoT, *tigT_bgoT_supp;

class SortCode {

	public :
		SortCode(){;} 
		void SortData(const char*, const char*, const char*, const char*, const char*);
		void Initialise();
};
#endif

void SortCode::Initialise() {

  printf("Start initializations\n");
  printf("Creating list\n");
  emmaList = new TList;
  icList = new TList;
  tigList = new TList;
  addList = new TList;
  tigtigList = new TList;
  tigbgoList = new TList;
  tigemmaList = new TList;
  addemmaList = new TList;
  timegatedList = new TList;
  massgatedList = new TList;
  ggatedList = new TList;
  ggatedicList = new TList;
  ggatedICSiList = new TList;
  printf("Creating histograms\n");

  xPos = new TH1F("xPos", "PGAC_X_Position", 160, -80, 80);
  emmaList->Add(xPos);
  yPos = new TH1F("yPos", "PGAC_Y_Position", 60, -30, 30);
  emmaList->Add(yPos);
  xPosT = new TH2F("xPosT", "PGAC_X_Position_V_Time", 4096, 0, 4096, 160, -80, 80);
  emmaList->Add(xPosT);
  yPosT = new TH2F("yPosT", "PGAC_Y_Position_V_Time", 4096, 0, 4096, 60, -30, 30);
  emmaList->Add(yPosT);
  pgac = new TH2F("pgac", "PGAC_Hit_Pattern", 160, -80, 80, 60, -30, 30);
  emmaList->Add(pgac);
  siE = new TH1F("siE", "EMMA_Focal_Plane_Silicon_Energy", 16384, 0, 16384);
  emmaList->Add(siE);
  siET = new TH2F("siET", "EMMA_Focal_Plane_Silicon_Energy_V_Time", 4096,0,4096, 16384, 0, 16384);
  emmaList->Add(siET);
  icSumVSi = new TH2F("icSumVSi","Si_Energy_VS_IC_SUM",4096,0,4096,16384,0,16384);
  emmaList->Add(icSumVSi);
  ic0VSi = new TH2F("ic0VSi","Si_Energy_VS_IC0",4096,0,4096,16384,0,16384);
  emmaList->Add(ic0VSi);
  ic1VSi = new TH2F("ic1VSi","Si_Energy_VS_IC1",4096,0,4096,16384,0,16384);
  emmaList->Add(ic1VSi);
  ic2VSi = new TH2F("ic2VSi","Si_Energy_VS_IC2",4096,0,4096,16384,0,16384);
  emmaList->Add(ic2VSi);
  ic3VSi = new TH2F("ic3VSi","Si_Energy_VS_IC3",4096,0,4096,16384,0,16384);
  emmaList->Add(ic3VSi);
  for (int i = 0; i < 5; i++) {
    char hname[64];
    sprintf(hname, "EMMA_IC-Segment_%i", i + 1);
    icE[i] = new TH1F(hname, hname, 4096, 0, 4096);
    icList->Add(icE[i]);
  }
  icN = new TH2F("icN", "IC_Energy_V_Segment_Number", 5, 0, 5, 4096, 0, 4096);
  icList->Add(icN);
  icSum = new TH1F("icSum", "IC_Summed_Energy", 16384, 0, 16384);
  icList->Add(icSum);
  iC0V1= new TH2F("iC0V1", "IC0_Energy_V_IC1_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC0V1);
  iC0V2= new TH2F("iC0V2", "IC0_Energy_V_IC2_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC0V2);
  iC0V3= new TH2F("iC0V3", "IC0_Energy_V_IC3_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC0V3);
  iC1V2= new TH2F("iC1V2", "IC1_Energy_V_IC2_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC1V2);
  iC1V3= new TH2F("iC1V3", "IC1_Energy_V_IC3_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC1V3);
  iC2V3= new TH2F("iC2V3", "IC2_Energy_V_IC3_Energy", 4096, 0, 4096, 4096, 0, 4096);
  icList->Add(iC2V3);

  for (int i = 0; i < 2; i++) {
    char hname[64];
    char hname2[64];
    sprintf(hname, "SSB_%i", i + 1);
    sprintf(hname2, "SSB_%i_Versus_Time", i + 1);
    ssbE[i] = new TH1F(hname, hname, 20000, 0, 200);
    ssbET[i] = new TH2F(hname2, hname2, 8000, 0, 8000, 2000, 0, 200);
    emmaList->Add(ssbE[i]);
    emmaList->Add(ssbET[i]);
  }

  tigE = new TH1F("tigE", "Tigress_Energy", 8192, 0, 8192);
  tigList->Add(tigE);
  tigDop = new TH1F("tigDop", "Tigress_Doppler_Corrected_Energy", 8192, 0, 8192);
  tigList->Add(tigDop);
  // begin CRN
  tigE_tigDop = new TH2F("tigE_tigDop", "Tigress_Energy_V_Doppler_Corrected_Energy", 8192, 0, 8192, 8192, 0, 8192);
  tigList->Add(tigE_tigDop);
  //end CRN
  tigE_ANum = new TH2F("tigE_ANum", "Tigress_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigE_ANum);
  tigDop_ANum = new TH2F("tigDop_ANum", "Tigress_Doppler_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigDop_ANum);
  tigE_theta = new TH2F("tigE_theta", "Tigress_Energy_V_Theta", 180, 0, 180, 8192, 0, 8192);
  tigList->Add(tigE_theta);
  tigDop_theta = new TH2F("tigDop_theta", "Tigress_Doppler_Energy_V_Theta", 180, 0, 180, 8192, 0, 8192);
  tigList->Add(tigDop_theta);

  addE = new TH1F("addE", "Addback_Energy", 8192, 0, 8192);
  addList->Add(addE);
  addDop = new TH1F("addDop", "Addback_Doppler_Corrected_Energy", 8192, 0, 8192);
  addList->Add(addDop);
  addE_ANum = new TH2F("addE_ANum", "Addback_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  addList->Add(addE_ANum);
  addDop_ANum = new TH2F("addDop_ANum", "Addback_Doppler_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  addList->Add(addDop_ANum);
  addE_theta = new TH2F("addE_theta", "Addback_Energy_V_Theta", 180, 0, 180, 8192, 0, 8192);
  addList->Add(addE_theta);
  addDop_theta = new TH2F("addDop_theta", "Addback_Doppler_Energy_V_Theta", 180, 0, 180, 8192, 0, 8192);
  addList->Add(addDop_theta);
  
  //TIGRESS-TIGRESS
  addT_addT = new TH1F("addT_addT","Tigress-Tigress_time",4096,-2048,2048); 
  tigtigList->Add(addT_addT);
  addE_addE = new TH2F("addE_addE", "Addback_Gamma_Gamma", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE);
  
  tigemmatof = new TH1F("tigemmatof", "Tigress_Time - EMMA Time", 20000, -10000, 10000);
  tigemmaList->Add(tigemmatof);
  addemmatof = new TH1F("addemmatof", "Addback_Time - EMMA Time", 20000, -10000, 10000);
  tigemmaList->Add(addemmatof);
  tigE_tof = new TH2F("tigE_tof", "Tigress_Energy_vs_(Tigress_Time - EMMA Time)", 20000, -10000, 10000, 8192, 0, 8192);
  tigemmaList->Add(tigE_tof);
  addE_tof = new TH2F("addE_tof", "Addback_Energy_vs_(Addback_Time - EMMA Time)", 20000, -10000, 10000, 8192, 0, 8192);
  tigemmaList->Add(addE_tof);
  tigDop_tof = new TH2F("tigDop_tof", "Tigress_Doppler_Energy_vs_(Tigress_Time - EMMA Time)", 20000, -10000, 10000, 8192, 0, 8192);
  tigemmaList->Add(tigDop_tof);
  addDop_tof = new TH2F("addDop_tof", "Addback_Doppler_Energy_vs_(Addback_Time - EMMA Time)", 20000, -10000, 10000, 8192, 0, 8192);
  tigemmaList->Add(addDop_tof);
  tDE = new TH2F("tDE","tDE",4000,0,4000,400,-2000,2000);
  tigemmaList->Add(tDE);
  
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
  
  tigDop_ANum_if_emma = new TH2F("tigDop_ANum_if_emma", "Dop_Energy_vs_Det_Num_if_emma", 64, 0, 64, 8192, 0, 8192);
  tigemmaList->Add(tigDop_ANum_if_emma);
  
  tigE_tg = new TH1F("tigE_tg", "Time_Gated_Tigress_Energy", 8192, 0, 8192);
  timegatedList->Add(tigE_tg);
  addDop_tg = new TH1F("addDop_tg", "Time_Gated_Addback_Doppler_Energy", 8192, 0, 8192);
  timegatedList->Add(addDop_tg);
  addDop_tg_bg = new TH1F("addDop_tg_bg", "Time_Gated_Addback_Doppler_Energy_Background", 8192, 0, 8192);
  timegatedList->Add(addDop_tg_bg);
  addDop_tg_bgs = new TH1F("addDop_tg_bgs", "Time_Gated_Addback_Doppler_Energy_BackgroundSubtracted", 8192, 0, 8192);
  addDop_tg_bgs->Sumw2();
  timegatedList->Add(addDop_tg_bgs);
  addDopp_addDopp_tg = new TH2F("addDopp_addDopp_tg", "Time_Gated_Addback_Doppler_Energy_gg", 4096, 0, 4096, 4096, 0, 4096);
  timegatedList->Add(addDopp_addDopp_tg);
  tigE_xPos_tg = new TH2F("tigE_xPos_tg", "Time_Gated_Tigress_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(tigE_xPos_tg);
  tigDop_xPos_tg = new TH2F("tigDopp_xPos_tg", "Time_Gated_Tigress_Doppler_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(tigDop_xPos_tg);
  tigE_xPos_bg = new TH2F("tigE_xPos_bg", "Time_Random_Tigress_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(tigE_xPos_bg);
  tigDop_xPos_bg = new TH2F("tigDop_xPos_bg", "Time_Random_Tigress_Doppler_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(tigDop_xPos_bg);
  addE_xPos_tg = new TH2F("addE_xPos_tg", "Time_Gated_Addback_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(addE_xPos_tg);
  addDop_xPos_tg = new TH2F("addDopp_xPos_tg", "Time_Gated_Addback_Doppler_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(addDop_xPos_tg);
  addE_xPos_bg = new TH2F("addE_xPos_bg", "Time_Random_Addback_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(addE_xPos_bg);
  addDop_xPos_bg = new TH2F("addDop_xPos_bg", "Time_Random_Addback_Doppler_Energy_vs_PGAC_X_position", 160, -80, 80, 8192, 0, 8192);
  timegatedList->Add(addDop_xPos_bg);
  pgac_tg = new TH2F("pgac_tg", "PGAC_Hit_Pattern_Time_Gated", 160, -80, 80, 60, -30, 30);
  timegatedList->Add(pgac_tg);
  pgac_bg = new TH2F("pgac_bg", "PGAC_Hit_Pattern_Time_Gated", 160, -80, 80, 60, -30, 30);
  timegatedList->Add(pgac_bg);
  icN_tg = new TH2F("icN_tg", "IC_Energy_V_Segment_Number_Time_Gated", 5, 0, 5, 4096, 0, 4096);
  timegatedList->Add(icN_tg);
  icN_bg = new TH2F("icN_bg", "IC_Energy_V_Segment_Number_Time_Random", 5, 0, 5, 4096, 0, 4096);
  timegatedList->Add(icN_bg);

  for (int i = 0; i < 2; i++) {
    char hname[64];
    char hname1[64];
    char hname2[64];
    char hname3[64];
    char hname4[64];
    char hname5[64];
    char hname6[64];
    char hname7[64];
    char hname8[64];
    sprintf(hname, "TIGRESS_Energy_MassGate_%i", i);
    sprintf(hname1, "TIGRESS_Doppler_Energy_MassGate_%i", i);
    sprintf(hname2, "Addback_Energy_MassGate_%i", i);
    sprintf(hname3, "Addback_Doppler_Energy_MassGate_%i", i);
    sprintf(hname4, "Time_Random_TIGRESS_Energy_MassGate_%i", i);
    sprintf(hname5, "Time_Random_TIGRESS_Doppler_Energy_MassGate_%i", i);
    sprintf(hname6, "Time_Random_Addback_Energy_MassGate_%i", i);
    sprintf(hname7, "Time_Random_Addback_Doppler_Energy_MassGate_%i", i);
    sprintf(hname8, "Addback_Doppler_ToF_MassGate_%i", i);
    tigE_massG[i] = new TH1F(hname, hname, 8192, 0, 8192);
    massgatedList->Add(tigE_massG[i]);
    tigDop_massG[i] = new TH1F(hname1, hname1, 8192, 0, 8192);
    massgatedList->Add(tigDop_massG[i]);
    addE_massG[i] = new TH1F(hname2, hname2, 8192, 0, 8192);
    massgatedList->Add(addE_massG[i]);
    addDop_massG[i] = new TH1F(hname3, hname3, 8192, 0, 8192);
    massgatedList->Add(addDop_massG[i]);
    tigE_massG_bg[i] = new TH1F(hname4, hname4, 8192, 0, 8192);
    massgatedList->Add(tigE_massG_bg[i]);
    tigDop_massG_bg[i] = new TH1F(hname5, hname5, 8192, 0, 8192);
    massgatedList->Add(tigDop_massG_bg[i]);
    addE_massG_bg[i] = new TH1F(hname6, hname6, 8192, 0, 8192);
    massgatedList->Add(addE_massG_bg[i]);
    addDop_massG_bg[i] = new TH1F(hname7, hname7, 8192, 0, 8192);
    massgatedList->Add(addDop_massG_bg[i]);
    addDop_mtof[i] = new TH1F(hname8, hname8, 20000, -10000, 10000);
    massgatedList->Add(addDop_mtof[i]);
  }

  for (int i = 0; i < 4; i++) {
    char hname[64];
    char hname1[64];
    char hname2[64];
    char hname3[64];
    char hname4[64];
    char hname5[64];

    char iname[64];
    char iname1[64];
    char iname2[64];
    char iname3[64];
    char iname4[64];
    char iname5[64];
    char iname6[64];
    char iname7[64];
    char iname8[64];
    char iname9[64];
    char iname10[64];
    char iname11[64];

    sprintf(hname, "PGAC_Hit_Pattern_Gamma_Gate_%i", i);
    sprintf(hname1, "PGAC_Hit_Pattern_Gamma_Gate_%i_Time_Random", i);
    sprintf(hname2, "PGAC_Hit_Pattern_Addback_Gate_%i", i);
    sprintf(hname3, "PGAC_Hit_Pattern_Addback_Gate_%i_Time_Random", i);
    sprintf(hname4, "IC_Energy_V_Segment_Number_Gamma_Gate_%i", i);
    sprintf(hname5, "IC_Energy_V_Segment_Number_Gamma_Gate_%i_Time_Random", i);

    sprintf(iname, "IC0_V_IC1_Gamma_Gate_%i", i);
    sprintf(iname1, "IC0_V_IC2_Gamma_Gate_%i", i);
    sprintf(iname2, "IC0_V_IC3_Gamma_Gate_%i", i);
    sprintf(iname3, "IC1_V_IC2_Gamma_Gate_%i", i);
    sprintf(iname4, "IC1_V_IC3_Gamma_Gate_%i", i);
    sprintf(iname10, "IC2_V_IC3_Gamma_Gate_%i", i);
    sprintf(iname5, "IC0_V_IC1_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname6, "IC0_V_IC2_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname7, "IC0_V_IC3_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname8, "IC1_V_IC2_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname9, "IC1_V_IC3_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname11, "IC2_V_IC3_Gamma_Gate_%i_Time_Random", i);

    pgac_gGate[i] = new TH2F(hname, hname, 160, -80, 80, 60, -30, 30);
    ggatedList->Add(pgac_gGate[i]);
    pgac_gGate_bg[i] = new TH2F(hname1, hname1, 160, -80, 80, 60, -30, 30);
    ggatedList->Add(pgac_gGate_bg[i]);
    pgac_aGate[i] = new TH2F(hname2, hname2, 160, -80, 80, 60, -30, 30);
    ggatedList->Add(pgac_aGate[i]);
    pgac_aGate_bg[i] = new TH2F(hname3, hname3, 160, -80, 80, 60, -30, 30);
    ggatedList->Add(pgac_aGate_bg[i]);
    icN_aGate_tg[i] = new TH2F(hname4, hname4, 5, 0, 5, 4096, 0, 4096);
    ggatedList->Add(icN_aGate_tg[i]);
    icN_aGate_bg[i] = new TH2F(hname5, hname5, 5, 0, 5, 4096, 0, 4096);
    ggatedList->Add(icN_aGate_bg[i]);

    iC0V1_aGate_tg[i] = new TH2F(iname, iname, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V1_aGate_tg[i]);
    iC0V2_aGate_tg[i] = new TH2F(iname1, iname1, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V2_aGate_tg[i]);
    iC0V3_aGate_tg[i] = new TH2F(iname2, iname2, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V3_aGate_tg[i]);
    iC1V2_aGate_tg[i] = new TH2F(iname3, iname3, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC1V2_aGate_tg[i]);
    iC1V3_aGate_tg[i] = new TH2F(iname4, iname4, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC1V3_aGate_tg[i]);
    iC2V3_aGate_tg[i] = new TH2F(iname10, iname10, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC2V3_aGate_tg[i]);

    iC0V1_aGate_bg[i] = new TH2F(iname5, iname5, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V1_aGate_bg[i]);
    iC0V2_aGate_bg[i] = new TH2F(iname6, iname6, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V2_aGate_bg[i]);
    iC0V3_aGate_bg[i] = new TH2F(iname7, iname7, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC0V3_aGate_bg[i]);
    iC1V2_aGate_bg[i] = new TH2F(iname8, iname8, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC1V2_aGate_bg[i]);
    iC1V3_aGate_bg[i] = new TH2F(iname9, iname9, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC1V3_aGate_bg[i]);
    iC2V3_aGate_bg[i] = new TH2F(iname11, iname11, 4096, 0, 4096, 4096, 0, 4096);
    ggatedicList->Add(iC2V3_aGate_bg[i]);
  }

  for (int i = 0; i < 4; i++) {
    char hname[64];
    char hname1[64];
    char hname2[64];
    char hname3[64];
    char hname4[64];

    char iname[64];
    char iname1[64];
    char iname2[64];
    char iname3[64];
    char iname4[64];

    sprintf(hname, "Si_Energy_VS_IC_SUM_Gamma_Gate_%i", i);
    sprintf(hname1, "Si_Energy_VS_IC0_Gamma_Gate_%i", i);
    sprintf(hname2, "Si_Energy_VS_IC1_Gamma_Gate_%i", i);
    sprintf(hname3, "Si_Energy_VS_IC2_Gamma_Gate_%i", i);
    sprintf(hname4, "Si_Energy_VS_IC3_Gamma_Gate_%i", i);

    sprintf(iname, "Si_Energy_VS_IC_SUM_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname1, "Si_Energy_VS_IC0_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname2, "Si_Energy_VS_IC1_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname3, "Si_Energy_VS_IC2_Gamma_Gate_%i_Time_Random", i);
    sprintf(iname4, "Si_Energy_VS_IC3_Gamma_Gate_%i_Time_Random", i);

    icSumVSi_aGate_tg[i] = new TH2F(hname, hname, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(icSumVSi_aGate_tg[i]);
    ic0VSi_aGate_tg[i] = new TH2F(hname1, hname1, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic0VSi_aGate_tg[i]);
    ic1VSi_aGate_tg[i] = new TH2F(hname2, hname2, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic1VSi_aGate_tg[i]);
    ic2VSi_aGate_tg[i] = new TH2F(hname3, hname3, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic2VSi_aGate_tg[i]);
    ic3VSi_aGate_tg[i] = new TH2F(hname4, hname4, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic3VSi_aGate_tg[i]);

    icSumVSi_aGate_bg[i] = new TH2F(iname, iname, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(icSumVSi_aGate_bg[i]);
    ic0VSi_aGate_bg[i] = new TH2F(iname1, iname1, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic0VSi_aGate_bg[i]);
    ic1VSi_aGate_bg[i] = new TH2F(iname2, iname2, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic1VSi_aGate_bg[i]);
    ic2VSi_aGate_bg[i] = new TH2F(iname3, iname3, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic2VSi_aGate_bg[i]);
    ic3VSi_aGate_bg[i] = new TH2F(iname4, iname4, 4096,0,4096,16384,0,16384);
    ggatedICSiList->Add(ic3VSi_aGate_bg[i]);
 } 
}
