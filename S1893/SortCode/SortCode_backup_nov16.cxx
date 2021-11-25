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
#include "TH3.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTrific.h"
#include "TS3Hit.h"
#include "TReaction.h"
#include "TSRIM.h"
#include "TTigress.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TSharc.h"
#include "TParserLibrary.h"
#include "TEnv.h"
#include "TCanvas.h"
#include "TApplication.h"
//#include "TSharcAnalysis.h"

using namespace std;




TList *tigList, *tigtigList, *tgList, *bgList, *PIDgatedList, *sharcList, *trifList, *trifHitmapList, *sharcTigList, *sharcProjList, *trifBraggList; 

TH1F *ssb_anode_dt1,*ssb_anode_dt2;

//Raw TIGRESS
TH1F *tigE, *addE, *addDoppRaw;
TH2F *tigE_ANum,*addE_ANum, *num_addr, *tigE_theta; 

//TIGRESS-TIGRESS
TH1F *addT_addT; 
TH2F *addE_addE, *addE_addE_tg;


//ToF Spectra
TH1F *excE_tg, *addDopp_tg, *excE_bg, *addDopp_bg, *addemmatof, *s3emmatof, *ssbICtof, *ssbSItof, *tigICtof, *tigAnodetof, *ssbEMTtof;
TH2F *addE_tof, *excE_theta_tg, *addDopp_ANum_tg, *addE_addE_tofg, *addDopp_exc_tg, *addDopp_addDopp_tg, *excE_theta_bg, *addDopp_ANum_bg, *addDopp_exc_bg, *hitmap_time_gated, *icSumVSi_gated_350, *addemmatofdet; 

//raw SHARC
TH1F *sharcE, *sharc_excE;
TH2F *sharcE_DetNum, *sharcE_front, *sharcE_back, *sharcE_phi, *sharcE_theta, *sharc_excE_theta, *sharc_dE_E, *sharc_XY_projection, *sharc_XZ_projection, *sharc_YZ_projection, *sharc_XY_projection_neg, *sharc_XZ_projection_neg, *sharc_YZ_projection_neg, *sharc_XY_projection_pos, *sharc_XZ_projection_pos, *sharc_YZ_projection_pos, *sharc_pad_vs_dssd;
TH3F *sharc_3D_hitmap;

TH2F *sharc_XY_array[17], *sharc_XZ_array[17], *sharc_YZ_array[17];

//raw TRIFIC
TH1F *trifE, *trif_XYStatus, *trifChan;
TH2F *trifE_chan, *trif_hitmap, *trif_BAD_hitmap, *trif_GOOD_hitmap, *trif_proj_YZ, *trif_proj_XZ, *trif_single_hitmap, *trif_braggE_detNum, *trifdEdx_Range, *trifdEdx_Range_interp_simple, *trifdEdx_Range_interp_complex;
TH3F *trif_3D_hitmap;

TH2F *trif_hitmap_array[33];
TGraph *trif_bragg_array[100];

//SHARC-TIGRESS
TH1F *tigSharcT, *tigE_sharc_gated;
TH2F *tigSharcTE, *sharcEtigE;




//TCuts
TCutG *massGate, *siicGate;

TCutG *na22_cut,*ne22_cut, *mg22_cut;
TCutG *na22_cut_ic,*ne22_cut_ic;

TCutG *pCut, *dCut, *tCut;

class SortCode {

	public :

		SortCode(){;} 
		void SortData(const char*, const char*, const char*, const char*);
		void Initialise();
};
#endif

void SortCode::Initialise() {

  printf("Start initializations\n");
  printf("Creating list\n");

  tigList = new TList;
  tigtigList = new TList;
  tgList = new TList;
  bgList = new TList;
  PIDgatedList = new TList;
  sharcList = new TList;
  trifList = new TList;
  trifHitmapList = new TList;
  sharcTigList = new TList;
  sharcProjList = new TList;
  trifBraggList = new TList;

  printf("Creating histograms\n");


  //Raw TIGRESS Spectra
  tigE = new TH1F("tigE", "Tigress_Energy", 8192, 0, 8192);
  tigList->Add(tigE);
  tigE_ANum = new TH2F("tigE_ANum", "Tigress_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigE_ANum);
  addE = new TH1F("addE", "Addback_Energy", 8192, 0, 8192);
  tigList->Add(addE);
  addDoppRaw = new TH1F("addDoppRaw", "Addback_Doppler_Raw", 8192, 0, 8192);
  tigList->Add(addDoppRaw);
  addE_ANum = new TH2F("addE_ANum", "Addback_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(addE_ANum);
  num_addr = new TH2F("num_addr", "TIGRESS_position_vs_address", 64, 0, 64, 16384, 0, 16384);
  tigList->Add(num_addr);
	tigE_theta = new TH2F("tigE_theta","TIGRESS energy vs lab angle",1024,0,2048,180,0,180);
	tigList->Add(tigE_theta);


 

  //TIGRESS-TIGRESS
  addT_addT = new TH1F("addT_addT","Tigress-Tigress_time",4096,-2048,2048); 
  tigtigList->Add(addT_addT);
  addE_addE = new TH2F("addE_addE", "Addback_Gamma_Gamma", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE);
  addE_addE_tg = new TH2F("addE_addE_tg", "Addback_Gamma_Gamma_Time_Gated", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE_tg);



 

  //Time Random Gated Spectra
  excE_bg = new TH1F("excE_bg","Excitation_Energy_MeV_Random_Time_Gated",2800,-2.,12.); 
  bgList->Add(excE_bg);
  excE_theta_bg = new TH2F("excE_theta_bg","Excitation_Energy_Theta_Random_Time_Gated",180,0,180,2800,-2.,12.); 
  bgList->Add(excE_theta_bg);
  addDopp_bg = new TH1F("addDopp_bg", "Addback_Doppler_Energy_Random_Time_Gated", 8192, 0, 8192);
  bgList->Add(addDopp_bg);
  addDopp_ANum_bg = new TH2F("addDopp_ANum_bg", "Addback_Doppler_Energy_V_ArrayNumber_Random_Time_Gated", 64, 0, 64, 8192, 0, 8192);
  bgList->Add(addDopp_ANum_bg);
  addDopp_exc_bg = new TH2F("addDopp_Exc_bg", "Addback_Doppler_Energy_V_Excitation_Energy_Random_Time_Gated", 8192, 0, 8192, 2800, -2., 12.);
  bgList->Add(addDopp_exc_bg);

//raw SHARC

  sharcE = new TH1F("sharcE","Sharc Energy; Energy (keV);Counts",8192,0,8192);
  sharcList->Add(sharcE);
  sharc_excE = new TH1F("sharc_excE","Sharc Excitation Energy;Energy (MeV);Counts",1000,0,20);
  sharcList->Add(sharc_excE);

  sharc_excE_theta = new TH2F("sharc_excE_theta","Sharc Excitation Energy vs Theta;#Theta (deg);Energy (MeV)",720,0,180,1000,0,20);
  sharcList->Add(sharc_excE_theta);
  sharcE_front = new TH2F("sharcE_front","Sharc Energy vs Front Strip;Strip;Energy (keV)",1700,0,1700,8192,0,8192);
  sharcList->Add(sharcE_front);
  sharcE_back = new TH2F("sharcE_back","Sharc Energy vs Back Strip;Strip;Energy (keV)",1700,0,1700,8192,0,8192);
  sharcList->Add(sharcE_back);
  sharcE_phi = new TH2F("sharcE_phi","Sharc Energy vs phi;#phi (deg);Energy (keV)",360,-180,180,8192,0,8192);
  sharcList->Add(sharcE_phi);
  sharcE_theta = new TH2F("sharcE_theta","Sharc Energy vs theta;#theta (deg);Energy (keV)",360,0,360,8192,0,8192);
  sharcList->Add(sharcE_theta);
  sharcE_DetNum = new TH2F("sharcE_DetNum","Sharc Energy vs Array Number;Det Number;Energy (keV)",20,0,20,8192,0,8192);
  sharcList->Add(sharcE_DetNum);
  sharc_dE_E = new TH2F("sharc_dE_E","Sharc dE vs E; Energy (keV); dE (keV)",3000,0,30000,1000,0,10000);
  sharcList->Add(sharc_dE_E);
  sharc_XY_projection = new TH2F("sharc_XY_projection","Sharc total XY Projection;X;Y",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XY_projection);
  sharc_XZ_projection = new TH2F("sharc_XZ_projection","Sharc total XZ Projection;X;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XZ_projection);
  sharc_YZ_projection = new TH2F("sharc_YZ_projection","Sharc total YZ Projection;Y;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_YZ_projection);

  sharc_XY_projection_neg = new TH2F("sharc_XY_projection_neg","Sharc XY Projection: Z < 0;X;Y",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XY_projection_neg);
  sharc_XZ_projection_neg = new TH2F("sharc_XZ_projection_neg","Sharc XZ Projection: Y < ;X;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XZ_projection_neg);
  sharc_YZ_projection_neg = new TH2F("sharc_YZ_projection_neg","Sharc YZ Projection: X < 0;Y;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_YZ_projection_neg);
  sharc_XY_projection_pos = new TH2F("sharc_XY_projection_pos","Sharc XY Projection: Z > 0;X;Y",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XY_projection_pos);
  sharc_XZ_projection_pos = new TH2F("sharc_XZ_projection_pos","Sharc XZ Projection: Y > 0;X;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XZ_projection_pos);
  sharc_YZ_projection_pos = new TH2F("sharc_YZ_projection_pos","Sharc YZ Projection: X > 0;Y;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_YZ_projection_pos);
  sharc_pad_vs_dssd = new TH2F("sharc_pad_vs_dssd","Sharc Pad # vs DSSD # ;DSSD # ;Pad #",20,0,20,20,0,20);
  sharcList->Add(sharc_pad_vs_dssd);


  sharc_3D_hitmap = new TH3F("sharc_3D_hitmap","Sharc 3D Hitmap;X;Y;Z",100,-100,100,100,-100,100,100,-100,100);
  sharcList->Add(sharc_3D_hitmap);
  
  //sharc projection maps for individual detectors
  for (int i = 0; i < 17; i++){
	sharc_XY_array[i] = new TH2F(Form("sharc_XY_det%i",i),Form("SHARC XY projection det %i;X;Y",i),100,-100,100,100,-100,100);
	sharcProjList->Add(sharc_XY_array[i]);
	sharc_XZ_array[i] = new TH2F(Form("sharc_XZ_det%i",i),Form("SHARC XZ projection det %i;X;Z",i),100,-100,100,100,-100,100);
	sharcProjList->Add(sharc_XZ_array[i]);
	sharc_YZ_array[i] = new TH2F(Form("sharc_YZ_det%i",i),Form("SHARC YZ projection det %i;Y;Z",i),100,-100,100,100,-100,100);
	sharcProjList->Add(sharc_YZ_array[i]);
  }
  

  

//raw TRIFIC
  trifE = new TH1F("trifE","Trific Energy;Energy (arb);Counts",1000,0,1000);
  trifList->Add(trifE);
  trifChan = new TH1F("trifChan","Trific Detector Number; Det. Number;Counts",35,-5,30);
  trifList->Add(trifChan);
  

  trifE_chan = new TH2F("trifE_chan","Trific Energy vs Grid Number;Grid Number;Energy (arb)",25,0,25,1000,0,1000);
  trifList->Add(trifE_chan); 
  trif_hitmap = new TH2F("trif_hitmap","Trific XY hitmap;X (mm);Y (mm);",200,-100,100,200,-100,100);
  trifList->Add(trif_hitmap);
  trif_single_hitmap = new TH2F("trif_single_hitmap","Trific XY Single Hitmap;X (mm);Y (mm)",200,-100,100,200,-100,100);
  trifList->Add(trif_single_hitmap);
  trif_BAD_hitmap = new TH2F("trif_BAD_hitmap","Trific XY hitmap (Bad events);X (mm);Y (mm);",200,-1000,1000,200,-1000,1000);
  trifList->Add(trif_BAD_hitmap);
  trif_GOOD_hitmap = new TH2F("trif_GOOD_hitmap","Trific XY hitmap (Good events);X (mm);Y (mm);",200,-100,100,200,-100,100);
  trifList->Add(trif_GOOD_hitmap);
  trif_proj_YZ = new TH2F("trif_proj_YZ","Trific YZ hitmap; Z (mm); Y(mm)",500,600,1100,200,-100,100);
  trifList->Add(trif_proj_YZ);
  trif_proj_XZ = new TH2F("trif_proj_XZ","Trific XZ hitmap; Z (mm); X(mm)",500,600,1100,200,-100,100);
  trifList->Add(trif_proj_XZ);
  trif_braggE_detNum = new TH2F("trif_braggE_detNum","Trific Bragg peak energy vs grid number;Grid Number;Bragg Peak Energy (arb.)",25,0,25,1000,0,1000);
  trifList->Add(trif_braggE_detNum);
  trifdEdx_Range = new TH2F("trifdEdx_Range","Trific dE/dx peak vs Range;Range (Grid Number);dE/dx (arb.)",35,0,35,1000,0,1000);
  trifList->Add(trifdEdx_Range);
  trifdEdx_Range_interp_simple = new TH2F("trifdEdx_Range_interp_simple","Trific dE/dx peak vs Range (simple interpolated);Range (cm);dE/dx (arb.)",70,0,70,300,0,300);
  trifList->Add(trifdEdx_Range_interp_simple);
  trifdEdx_Range_interp_complex = new TH2F("trifdEdx_Range_interp_complex","Trific dE/dx peak vs Range (complex interpolated);Range (cm);dE/dx (arb.)",70,0,70,300,0,300);
  trifList->Add(trifdEdx_Range_interp_complex);


  trif_XYStatus = new TH1F("trif_XYStatus","Trific XY Reconstruction status",5,0,5);
  TAxis *ax = trif_XYStatus->GetXaxis();
  ax->SetBinLabel(2,"Good XY Reconstruction");
  ax->SetBinLabel(4,"Bad XY Reconstruction");
  trifList->Add(trif_XYStatus);

  trif_3D_hitmap = new TH3F("trif_3D_hitmap","Trific XYZ Hitmap;X (mm); Y (mm); Z(mm)",200,-100,100,200,-100,100,500,600,1100);
  trifList->Add(trif_3D_hitmap);


// TRIFIC Hitmaps
  for (int i = 0; i < 33; i++){ //we do 0-32 even though trific only has grids 1-24 b/c some unused channels could be still cables up and in the odb.
     trif_hitmap_array[i] = new TH2F(Form("trif_hitmap_%i",i),Form("Trific XY hitmap: Grid %i;X (mm);Y (ymm)",i),200,-100,100,200,-100,100);
     trifHitmapList->Add(trif_hitmap_array[i]);
  }
  /*for (int i = 0; i < 100; i++){ //this will create 100 trific bragg graphs to use for debugging.
	trif_bragg_array[i] = new TGraph(Form("trif_bragg_%i",i),Form("Trific Bragg Graph %i; Grid Number;dE/dx (arb.)",30,0,30,1000,0,1000);
	trifBraggList->Add(trif_bragg_array[i]);
  }*/
  

//SHARC-TIGRESS
  tigSharcT = new TH1F("tigT_sharcT","TIGRESS-SHARC timing;ttig - tsharc;Counts",8192,-4096,4096);
  sharcTigList->Add(tigSharcT);
  tigSharcTE = new TH2F("tigSharcTE","TIGRESS-SHARC timing vs TIGRESS energy;ttig - tsharc;TIGRESS Energy (keV);",8192,-4096,4096,8192,0,8192);
  sharcTigList->Add(tigSharcTE);
  tigE_sharc_gated = new TH1F("tigE_sharc_gated","TIGRESS Energy in coinc. with SHARC hit;TIG Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_gated);
  sharcEtigE = new TH2F("sharcEtigE","TIGRESS energy vs SHARC energy;SHARC Energy (keV);TIGRESS Energy (keV);",8192,0,8192,8192,0,8192);
  sharcTigList->Add(sharcEtigE);
  

//Loading TCuts
  printf("\n\nLoading cut files:\n");
  TFile *f1 = new TFile("/tig/pterodon_data3/S1893/SortCode/Cuts/pCut.root"); //protons
  pCut = (TCutG*)f1->FindObjectAny("pCut");
  f1->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/dCut.root"); //deuterons
  dCut = (TCutG*)f1->FindObjectAny("dCut");
  f1->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/tCut.root"); //tritons
  dCut = (TCutG*)f1->FindObjectAny("tCut");
  f1->Close();

  

}
