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
#include "TMultiGraph.h"
//#include "TSharcAnalysis.h"

using namespace std;




TList *tigList, *tigtigList, *tgList, *bgList, *PIDgatedList, *sharcList, *trifList, *trifHitmapList, *sharcTigList, *trifTigList, *sharcArrayList, *trifBraggList, *sharcTrificList, *trifPIDGatedList; 

TH1F *ssb_anode_dt1,*ssb_anode_dt2;

//Raw TIGRESS
TH1F *tigE, *addE, *addDoppRaw, *tigDopp;
TH2F *tigE_ANum,*addE_ANum, *num_addr, *tigE_theta, *tigDopp_ANum; 

//TIGRESS-TIGRESS
TH1F *addT_addT; 
TH2F *addE_addE, *addE_addE_tg;


//ToF Spectra
TH1F *excE_tg, *addDopp_tg, *excE_bg, *addDopp_bg, *addemmatof, *s3emmatof, *ssbICtof, *ssbSItof, *tigICtof, *tigAnodetof, *ssbEMTtof;
TH2F *addE_tof, *excE_theta_tg, *addDopp_ANum_tg, *addE_addE_tofg, *addDopp_exc_tg, *addDopp_addDopp_tg, *excE_theta_bg, *addDopp_ANum_bg, *addDopp_exc_bg, *hitmap_time_gated, *icSumVSi_gated_350, *addemmatofdet; 

//raw SHARC
TH1F *sharcE, *sharc_excE;
TH2F *sharcE_FB;
TH1F *sharcT_FB;
TH2F *sharcE_DetNum, *sharcE_front, *sharcE_back, *sharcE_phi, *sharcE_theta, *sharc_excE_theta, *sharc_dE_E, *sharc_XY_projection, *sharc_XZ_projection, *sharc_YZ_projection, *sharc_XY_projection_neg, *sharc_XZ_projection_neg, *sharc_YZ_projection_neg, *sharc_XY_projection_pos, *sharc_XZ_projection_pos, *sharc_YZ_projection_pos, *sharc_pad_vs_dssd;
TH3F *sharc_3D_hitmap;
TH2F *sharc_ffE, *sharc_bfE, *sharcEB_theta;

TH2F *sharcDebug_dE_E, *sharcDebug_dE_PadE, *sharcDebug_E_PadE, *sharcDebug_dEPadE_E;


TH2F *sharc_XY_array[17], *sharc_XZ_array[17], *sharc_YZ_array[17], *sharcE_theta_array[17], *sharc_dE_E_array[9];

//raw TRIFIC
TH1F *trifE, *trif_XYStatus, *trifChan;
TH2F *trifEdE;
TH2F *trifE_chan, *trif_hitmap, *trif_BAD_hitmap, *trif_GOOD_hitmap, *trif_proj_YZ, *trif_proj_XZ, *trif_single_hitmap, *trif_braggE_detNum, *trifdEdx_Range, *trifdEdx_Range_interp_simple, *trifdEdx_Range_interp_complex, *trifE_chan_cut, *trif_PID;
TH3F *trif_3D_hitmap;

TH2F *trifEdE_93Sr_accepted, *tridEdE_93Rb_accepted;


TH2F *trif_hitmap_array[33];
TGraph *trif_bragg_array[100];

TMultiGraph *trifBraggMulti;

//SHARC-TIGRESS
TH1F *tigSharcT, *tigE_sharc_gated, *sharcE_tig_gated, *tigE_sharc_p_gated, *tigE_sharc_d_gated, *tigE_sharc_t_gated, *tigE_sharc_a_gated, *tigE_sharc_c_gated;
TH2F *tigSharcTE, *sharcEtigE, *tigSharcT_front, *tigSharcT_back, *tigSharcT_phi;
TH1F *tigDop_gated;
TH2F *tigDop_gated_scan;
TH1F *tig_sharc_cut_dopp;
TH1F *tigDop_gated_bkgn_low, *tigDop_gated_bkgn_high;

//SHARC-TRIFIC
TH1F *trificSharcT, *trificSharcT_GOOD;

//TRIFIC-TIG
TH1F *trif_cut1_TigE, *trif_cut2_TigE, *TTTime, *trif_cut_93SrPID_tigE, *trif_cut_93RbPID_tigE, *trif_cut_93SrPID_DoppE, *trif_cut_93RbPID_DoppE;
TH2F *trif_tig_gated_scan;

//TCuts
TCutG *massGate, *siicGate;

TCutG *na22_cut,*ne22_cut, *mg22_cut;
TCutG *na22_cut_ic,*ne22_cut_ic;

TCutG *pCut, *dCut, *tCut, *aCut, *cCut; //proton, deuteron, triton, alpha, carbon cuts on dE vs E curve (general)
TCutG *pCut_array[9], *dCut_array[9], *tCut_array[9], *aCut_array[9], *cCut_array[9];//individual detector cuts
TCutG *cutg;

//trific pid cuts
TCutG *trifCut1, *trifCut2, *trifCut3;
TCutG *sr93PIDCut, *rb93PIDCut;



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
  sharcTrificList = new TList;
  sharcArrayList = new TList;
  trifBraggList = new TList;
  trifPIDGatedList = new TList;
  trifTigList = new TList;

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
	
  tigDopp = new TH1F("tigDopp", "Tigress_Doppler_Energy", 8192, 0, 8192);
  tigList->Add(tigDopp);
  tigDopp_ANum = new TH2F("tigDopp_ANum", "Tigress_Doppler_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigDopp_ANum);	 


 

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
  sharc_excE = new TH1F("sharc_excE","Sharc Excitation Energy;Energy (MeV);Counts",2800,-2,12);
  sharcList->Add(sharc_excE);

  sharc_excE_theta = new TH2F("sharc_excE_theta","Sharc Excitation Energy vs Theta;#Theta (deg);Energy (MeV)",720,0,180,1400,-2,12);
  sharcList->Add(sharc_excE_theta);

  sharcE_FB = new TH2F("sharcE_FB","Sharc front vs back energy",1024,0,8192,1024,0,8192);
  sharcList->Add(sharcE_FB);
  sharcT_FB = new TH1F("sharcT_FB","Sharc front back time difference",1000,-2000,2000);
  sharcList->Add(sharcT_FB);
  sharcE_front = new TH2F("sharcE_front","Sharc Energy vs Front Strip;Strip;Energy (keV)",1700,0,1700,8192,0,8192);
  sharcList->Add(sharcE_front);
  sharcE_back = new TH2F("sharcE_back","Sharc Energy vs Back Strip;Strip;Energy (keV)",1700,0,1700,8192,0,8192);
  sharcList->Add(sharcE_back);
  sharcE_phi = new TH2F("sharcE_phi","Sharc Energy vs phi;#phi (deg);Energy (keV)",360,-180,180,8192,0,8192);
  sharcList->Add(sharcE_phi);
  sharcE_theta = new TH2F("sharcE_theta","Sharc Energy vs theta;#theta (deg);Energy (keV)",180,0,180,16384,0,16384);
  sharcList->Add(sharcE_theta);
  sharcEB_theta =  new TH2F("sharcEB_theta","Sharc Energy vs theta;theta;Energy",180,0,180,16384,0,16384);
  sharcList->Add(sharcEB_theta);
  sharcE_DetNum = new TH2F("sharcE_DetNum","Sharc Energy vs Array Number;Det Number;Energy (keV)",20,0,20,16384,0,16384);
  sharcList->Add(sharcE_DetNum);
  sharc_dE_E = new TH2F("sharc_dE_E","Sharc dE vs E; Energy (keV); dE (arb)",199,100,2e4,3e2,0,3e4);
  sharcList->Add(sharc_dE_E);
  sharc_XY_projection = new TH2F("sharc_XY_projection","Sharc total XY Projection;X;Y",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XY_projection);
  sharc_XZ_projection = new TH2F("sharc_XZ_projection","Sharc total XZ Projection;X;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_XZ_projection);
  sharc_YZ_projection = new TH2F("sharc_YZ_projection","Sharc total YZ Projection;Y;Z",100,-100,100,100,-100,100);
  sharcList->Add(sharc_YZ_projection);
  sharc_ffE = new TH2F("sharc_ffE","Sharc Front Fragment Energy",900,0,900,8192,0,8192);
  sharc_bfE = new TH2F("sharc_bfE","Sharc Back Fragment Energy",900,0,900,8192,0,8192);
  sharcList->Add(sharc_ffE);
  sharcList->Add(sharc_bfE);
  
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
  sharc_pad_vs_dssd = new TH2F("sharc_pad_vs_dssd","Sharc Pad # vs DSSD # ;DSSD # ;Pad #",21,-1,20,21,-1,20);
  sharcList->Add(sharc_pad_vs_dssd);

  sharc_3D_hitmap = new TH3F("sharc_3D_hitmap","Sharc 3D Hitmap;X;Y;Z",100,-100,100,100,-100,100,100,-100,100);
  sharcList->Add(sharc_3D_hitmap);
  
  //sharc projection maps for individual detectors
  for (int i = 0; i < 17; i++){
	sharc_XY_array[i] = new TH2F(Form("sharc_XY_det%i",i),Form("SHARC XY projection det %i;X;Y",i),100,-100,100,100,-100,100);
	sharcArrayList->Add(sharc_XY_array[i]);
	sharc_XZ_array[i] = new TH2F(Form("sharc_XZ_det%i",i),Form("SHARC XZ projection det %i;X;Z",i),100,-100,100,100,-100,100);
	sharcArrayList->Add(sharc_XZ_array[i]);
	sharc_YZ_array[i] = new TH2F(Form("sharc_YZ_det%i",i),Form("SHARC YZ projection det %i;Y;Z",i),100,-100,100,100,-100,100);
	sharcArrayList->Add(sharc_YZ_array[i]);
	sharcE_theta_array[i] = new TH2F(Form("sharcE_theta_det%i",i),Form("SHARC Energy vs #theta det %i;#theta (deg);Energy (keV)",i),360,0,360,16384,0,16384);
	sharcArrayList->Add(sharcE_theta_array[i]);

  }

  for (int i = 0; i < 9; i++){
	sharc_dE_E_array[i] = new TH2F(Form("sharc_dE_E_det%i",i),Form("SHARC #DeltaE vs E det %i;E (arb);#DeltaE (arb)",i),199,100,2e4,3e2,0,3e4);
	sharcArrayList->Add(sharc_dE_E_array[i]);
  }


  sharcDebug_dE_E = new TH2F("sharcDebug_dE_E","SHARC Debugging dE vs E;E (arb);dE (arb)",199,100,2e4,3e2,0,3e4);
  sharcList->Add(sharcDebug_dE_E);
  sharcDebug_dE_PadE = new TH2F("sharcDebug_dE_PadE","SHARC Debugging dE vs Pad E;Pad E (arb);dE (arb)",199,100,2e4,3e2,0,3e4);
  sharcList->Add(sharcDebug_dE_PadE);
  sharcDebug_E_PadE = new TH2F("sharcDebug_E_PadE","SHARC Debugging E vs Pad E;Pad E (arb);E (arb)",199,100,2e4,3e2,0,3e4);
  sharcList->Add(sharcDebug_E_PadE);
  sharcDebug_dEPadE_E = new TH2F("sharcDebug_dEPadE_E","SHARC Debugging dE+padE vs total E;E (arb);dE+padE (arb)",199,100,2e4,3e2,0,3e4);
  sharcList->Add(sharcDebug_dEPadE_E);
  

  

//raw TRIFIC
  trifEdE = new TH2F("trifEdE","Trific Energy - dE;dE (arb);E (arb)",2048,0,32768,2048,0,32768);
  trifList->Add(trifEdE);
  trifE = new TH1F("trifE","Trific Energy;Energy (arb);Counts",10000,0,10000);
  trifList->Add(trifE);
  trifChan = new TH1F("trifChan","Trific Detector Number; Det. Number;Counts",35,-5,30);
  trifList->Add(trifChan);
  

  trifE_chan = new TH2F("trifE_chan","Trific Energy vs Grid Number;Grid Number;Energy (arb)",25,0,25,10000,0,10000);
  trifList->Add(trifE_chan); 
  trifE_chan_cut = new TH2F("trifE_chan_cut","Trific Energy vs Grid Number Gated;Grid Number;Energy (arb)",25,0,25,10000,0,10000);
  trifList->Add(trifE_chan_cut); 
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
  trifdEdx_Range = new TH2F("trifdEdx_Range","Trific dE/dx peak vs Range;Range (Grid Number);dE/dx (arb.)",35,0,35,10000,0,10000);
  trifList->Add(trifdEdx_Range);
  trifdEdx_Range_interp_simple = new TH2F("trifdEdx_Range_interp_simple","Trific dE/dx peak vs Range (simple interpolated);Range (cm);dE/dx (arb.)",200,0,200,3000,0,3000);
  trifList->Add(trifdEdx_Range_interp_simple);
  trifdEdx_Range_interp_complex = new TH2F("trifdEdx_Range_interp_complex","Trific dE/dx peak vs Range (complex interpolated);Range (cm);dE/dx (arb.)",200,0,200,3000,0,3000);
  trifList->Add(trifdEdx_Range_interp_complex);

  trifEdE_93Sr_accepted = new TH2F("trifEdE_93Sr_accepted","Accepted 93Sr: Trific Energy - dE;Energy (arb);Counts",2048,0,32768,2048,0,32768);
  trifList->Add(trifEdE_93Sr_accepted);
  tridEdE_93Rb_accepted = new TH2F("tridEdE_93Rb_accepted","Accepted 93Rb: Trific Energy - dE;Energy (arb);Counts",2048,0,32768,2048,0,32768);
  trifList->Add(tridEdE_93Rb_accepted);

  trif_PID = new TH2F("trif_PID","Trific PID plot;E (arb);dE (arb)",5000,0,20000,5000,0,20000);
  trifList->Add(trif_PID);



  trif_XYStatus = new TH1F("trif_XYStatus","Trific XY Reconstruction status",5,0,5);
  TAxis *ax = trif_XYStatus->GetXaxis();
  ax->SetBinLabel(2,"Good XY Reconstruction");
  ax->SetBinLabel(4,"Bad XY Reconstruction");
  trifList->Add(trif_XYStatus);

  trif_3D_hitmap = new TH3F("trif_3D_hitmap","Trific XYZ Hitmap;X (mm); Y (mm); Z(mm)",200,-100,100,200,-100,100,500,600,1100);
  trifList->Add(trif_3D_hitmap);


  trifBraggMulti = new TMultiGraph();
  trifList->Add(trifBraggMulti);


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
  //TH1F
  tigSharcT = new TH1F("tigT_sharcT","TIGRESS-SHARC timing;ttig - tsharc;Counts",8192,-4096,4096);
  sharcTigList->Add(tigSharcT);
  tigE_sharc_gated = new TH1F("tigE_sharc_gated","TIGRESS Energy in coinc. with SHARC hit;TIG Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_gated);
  sharcE_tig_gated = new TH1F("sharcE_tig_gated","Sharc Energy in coinc. with TIGRESS 60 keV; Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(sharcE_tig_gated);
  tig_sharc_cut_dopp = new TH1F("tig_sharc_cut_dopp","Tig gated on Sharccut Dopp cor",4000,0,4000);
  sharcTigList->Add(tig_sharc_cut_dopp);
  tigDop_gated = new TH1F("tigDop_gated","Tig Dop En gated on SharcEn and dT",8192,0,8192);
  sharcTigList->Add(tigDop_gated);
  tigE_sharc_p_gated = new TH1F("tigE_sharc_p_gated","Tig Eng gated on sharc protons;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_p_gated);
  tigE_sharc_d_gated = new TH1F("tigE_sharc_d_gated","Tig Eng gated on sharc deuterons;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_d_gated);
  tigE_sharc_t_gated = new TH1F("tigE_sharc_t_gated","Tig Eng gated on sharc tritons;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_t_gated);
  tigE_sharc_a_gated = new TH1F("tigE_sharc_a_gated","Tig Eng gated on sharc alphas;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_a_gated);
  tigE_sharc_c_gated = new TH1F("tigE_sharc_c_gated","Tig Eng gated on sharc carbons?;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigE_sharc_c_gated);

  tigDop_gated_bkgn_low = new TH1F("tigDop_gated_bkgn_low","Tig dopp eng background low;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigDop_gated_bkgn_low);
  tigDop_gated_bkgn_high = new TH1F("tigDop_gated_bkgn_high","Tig dopp eng background high;Energy (keV);Counts",8192,0,8192);
  sharcTigList->Add(tigDop_gated_bkgn_high);
  


  //TH2F
  tigSharcTE = new TH2F("tigSharcTE","TIGRESS-SHARC timing vs TIGRESS energy;ttig - tsharc;TIGRESS Energy (keV);",8192,-4096,4096,8192,0,8192);
  sharcTigList->Add(tigSharcTE);
  tigSharcT_front = new TH2F("tigSharcT_front","TIGRESS-SHARC timing vs Front Strip;Strip;ttig - tsharc",1700,0,1700,8192,-4096,4096);
  sharcTigList->Add(tigSharcT_front);
  tigSharcT_back = new TH2F("tigSharcT_back","TIGRESS-SHARC timing vs Back Strip;Strip;ttig - tsharc",1700,0,1700,8192,-4096,4096);
  sharcTigList->Add(tigSharcT_back);
  tigSharcT_phi = new TH2F("tigSharcT_phi","TIGRESS-SHARC timing vs phi;#phi (deg);ttig - tsharc",360,-180,180,8192,-4096,4096);
  sharcTigList->Add(tigSharcT_phi);
  sharcEtigE = new TH2F("sharcEtigE","TIGRESS energy vs SHARC energy;SHARC Energy (keV);TIGRESS Energy (keV);",8192,0,8192,8192,0,8192);
  sharcTigList->Add(sharcEtigE);
  tigDop_gated_scan = new TH2F("tigDop_gated_scan","Tig Dop En gated on SharcEn and dT vs Beta",20,0,20,4000,0,4000);
  sharcTigList->Add(tigDop_gated_scan);

  
  
  
   
	

//SHARC-TRIFIC
  trificSharcT = new TH1F("trificT_sharcT","TRIFIC-SHARC timing;ttrific - tsharc;Counts",8192,-4096,4096);
  sharcTrificList->Add(trificSharcT);
  trificSharcT_GOOD = new TH1F("trificT_sharcT_GOOD","TRIFIC-SHARC timing (Good events);ttrific - tsharc;Counts",8192,-4096,4096);
  sharcTrificList->Add(trificSharcT_GOOD);

//Loading TCuts
  printf("\n\nLoading cut files:\n");
 /* TFile *f1 = new TFile("/tig/pterodon_data3/S1893/SortCode/Cuts/pCut.root"); //protons
  pCut = (TCutG*)f1->FindObjectAny("pCut");
  pCut = (TCutG*)pCut->Clone();
  f1->Close();
  TFile *f2 = f2->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/dCut.root"); //deuterons
  dCut = (TCutG*)f2->FindObjectAny("dCut");
  dCut = (TCutG*)dCut->Clone();
  f2->Close();
  TFile *f3 = f3->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/tCut.root"); //tritons
  tCut = (TCutG*)f3->FindObjectAny("tCut");
  tCut = (TCutG*)tCut->Clone();
  f3->Close();
 //TFile *f4 = f4->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/aCut.root"); //alphas
  //aCut = (TCutG*)f4->FindObjectAny("aCut");
  //aCut = (TCutG*)aCut->Clone();
  //f4->Close();
  TFile *f5 = f5->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/cCut.root"); //carbon
  cCut = (TCutG*)f5->FindObjectAny("cCut");
  cCut = (TCutG*)cCut->Clone();
  f5->Close();*/


//TRIFIC-TIG
  trif_cut1_TigE = new TH1F("trif_cut1_TigE","TIGRESS E gated on TRIFIC Cut 1;Energy (keV);Counts",4000,0,4000);
  //trifPIDGatedList->Add(trif_cut1_TigE);
  trif_cut2_TigE = new TH1F("trif_cut2_TigE","TIGRESS E gated on TRIFIC Cut 2;Energy (keV);Counts",4000,0,4000);
  //trifPIDGatedList->Add(trif_cut2_TigE);
  TTTime = new TH1F("TTTime","Triffic-Tigress time difference",10000,-10000,10000);
  trif_tig_gated_scan = new TH2F("trif_tig_gated_scan","TigDop triffic gate vs beta",20,0,20,4000,0,4000);
  trifTigList->Add(TTTime);
  trifTigList->Add(trif_cut1_TigE);
  trifTigList->Add(trif_cut2_TigE);
  trifTigList->Add(trif_tig_gated_scan);

  trif_cut_93SrPID_tigE = new TH1F("trif_cut_93SrPID_tigE","TIGRESS E Gated on 93Sr PID;Energy (keV);Counts",8192,0,8192);
  trifTigList->Add(trif_cut_93SrPID_tigE);
  trif_cut_93RbPID_tigE = new TH1F("trif_cut_93RbPID_tigE","TIGRESS E Gated on 93Rb PID;Energy (keV);Counts",8192,0,8192);
  trifTigList->Add(trif_cut_93RbPID_tigE);
  trif_cut_93SrPID_DoppE = new TH1F("trif_cut_93SrPID_DoppE","TIGRESS Doppler E Gated on 93Sr PID;Energy (keV);Counts",8192,0,8192);
  trifTigList->Add(trif_cut_93SrPID_DoppE);
  trif_cut_93RbPID_DoppE = new TH1F("trif_cut_93RbPID_DoppE","TIGRESS Doppler E Gated on 93Rb PID;Energy (keV);Counts",8192,0,8192);
  trifTigList->Add(trif_cut_93RbPID_DoppE);


  
  for (int i = 5; i < 9; i++){
	TCutG *temp;
	TFile *f6 = f6->Open(Form("/tig/pterodon_data3/S1893/SortCode/Cuts/pCut_%i.root",i)); //protons
	temp = (TCutG*)f6->FindObjectAny("pCut");
	pCut_array[i] = (TCutG*)temp->Clone();
	f6->Close();
	TFile *f7 = f7->Open(Form("/tig/pterodon_data3/S1893/SortCode/Cuts/dCut_%i.root",i)); //deutrons
	dCut_array[i] = (TCutG*)f7->FindObjectAny("dCut");
	dCut_array[i] = (TCutG*)dCut_array[i]->Clone();
	f7->Close();
	TFile *f8 = f8->Open(Form("/tig/pterodon_data3/S1893/SortCode/Cuts/tCut_%i.root",i)); //tritons
	tCut_array[i] = (TCutG*)f8->FindObjectAny("tCut");
	tCut_array[i] = (TCutG*)tCut_array[i]->Clone();
	f8->Close();
	//TFile *f9 = f9->Open(Form("/tig/pterodon_data3/S1893/SortCode/Cuts/aCut_%i.root",i)); //alphas
	//aCut_array[i] = (TCutG*)f9->FindObjectAny("aCut");
	//aCut_array[i] = (TCutG*)aCut_array[i]->Clone();
	//f9->Close();*/
	TFile *f10 = f10->Open(Form("/tig/pterodon_data3/S1893/SortCode/Cuts/cCut_%i.root",i)); //carbons
	cCut_array[i] = (TCutG*)f10->FindObjectAny("cCut");
	cCut_array[i] = (TCutG*)cCut_array[i]->Clone();
	f10->Close();
  }

  
   /*TFile *trif1 = trif1->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/trifCut1.root"); //first cut
   trifCut1 = (TCutG*)trif1->FindObjectAny("trifCut");
   trifCut1 = (TCutG*)trifCut1->Clone();
   trif1->Close();

   TFile *trif2 = trif2->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/trifCut2.root"); //second cut
   trifCut2 = (TCutG*)trif2->FindObjectAny("trifCut");
   trifCut2 = (TCutG*)trifCut2->Clone();
   trif2->Close();

   TFile *trif3 = trif3->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/trifCut3.root"); //third cut
   trifCut3 = (TCutG*)trif3->FindObjectAny("trifCut");
   trifCut3 = (TCutG*)trifCut3->Clone();
   trif3->Close();
*/

   TFile *srPIDFile = srPIDFile->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/Sr93PIDCut.root");
   sr93PIDCut = (TCutG*)srPIDFile->FindObjectAny("Sr93Cut");
   sr93PIDCut = (TCutG*)sr93PIDCut->Clone();
   srPIDFile->Close();

   TFile *rbPIDFile = rbPIDFile->Open("/tig/pterodon_data3/S1893/SortCode/Cuts/Rb93PIDCut.root");
   rb93PIDCut = (TCutG*)rbPIDFile->FindObjectAny("Rb93Cut");
   rb93PIDCut = (TCutG*)rb93PIDCut->Clone();
   rbPIDFile->Close();


   cutg = new TCutG("pCut",10);
   cutg->SetVarX("SHARC #DeltaE vs E det 7");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,194.444,3059.62);
   cutg->SetPoint(1,1250,1746.27);
   cutg->SetPoint(2,3361.11,1113.91);
   cutg->SetPoint(3,7361.11,822.059);
   cutg->SetPoint(4,8916.67,1065.27);
   cutg->SetPoint(5,5361.11,1527.38);
   cutg->SetPoint(6,2861.11,2135.41);
   cutg->SetPoint(7,1277.78,2986.65);
   cutg->SetPoint(8,333.333,3862.22);
   cutg->SetPoint(9,194.444,3059.62);
   cutg->Draw("");
  

}
