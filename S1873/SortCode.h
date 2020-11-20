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
#include "TS3.h"
#include "TS3Hit.h"
#include "TReaction.h"
#include "TSRIM.h"
#include "TTigress.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TEmma.h"
#include "TParserLibrary.h"
#include "TEnv.h"
#include "Declarations.h"
using namespace std;

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

  tigList = new TList;
  s3List = new TList;
  tigtigList = new TList;
  tigs3List = new TList;
  emmaList = new TList;
  icList = new TList;
  tigemmaList = new TList;
  tgList = new TList;
  bgList = new TList;
  PIDgatedList = new TList;
  PIDgatedListS3 = new TList;

  printf("Creating histograms\n");


  //Raw TIGRESS Spectra
  tigE = new TH1F("tigE", "Tigress_Energy", 8192, 0, 8192);
  tigList->Add(tigE);
  tigE_ANum = new TH2F("tigE_ANum", "Tigress_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(tigE_ANum);
  addE = new TH1F("addE", "Addback_Energy", 8192, 0, 8192);
  tigList->Add(addE);
  addE_ANum = new TH2F("addE_ANum", "Addback_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigList->Add(addE_ANum);
  num_addr = new TH2F("num_addr", "TIGRESS_position_vs_address", 64, 0, 64, 16384, 0, 16384);
  tigList->Add(num_addr);

  //Raw S3 Spectra
  rings = new TH2F("rings","rings",24,0,24,2048,0,8192); 
  s3List->Add(rings);
  sectors = new TH2F("sectors","sectors",32,0,32,2048,0,8192); 
  s3List->Add(sectors);
  s3_E = new TH1F("S3_E","S3_E",8192,0,16384);
  s3List->Add(s3_E); 
  s3_E_theta = new TH2F("s3_E_theta","s3_E_theta",180,0,180, 8192, 0, 16384); 
  s3List->Add(s3_E_theta);
  hitmap = new TH2F("hitmap","Hitmap",200,-50,50,200,-50,50); 
  s3List->Add(hitmap);
  s3_rings_sectors = new TH2F("s3_rings_sectors","S3_sectorE_vs_ringE;Ring Energy;Sector Energy",2048,0,8192,2048,0,8192);
  s3List->Add(s3_rings_sectors);
  s3_rings_sectors_singles = new TH2F("s3_rings_sectors_singles","S3_sectorE_vs_ringE_Multiplicity_1;Ring Energy;Sector Energy",2048,0,8192,2048,0,8192);
  s3List->Add(s3_rings_sectors_singles);
  s3_rings_sectors_singlesT = new TH1F("s3_rings_sectors_singles_timing","S3_sectorT_vs_ringT_Multiplicity_1;Sector Time - Ring Time",2048,-1024,1024);
  s3List->Add(s3_rings_sectors_singlesT);
  s3_rings_sectors_singlesTvT = new TH2F("s3_rings_sectors_singles_timing_v_time","S3_sectorT_vs_ringT_Multiplicity_1_vs_Time;Time (min);Sector Time - Ring Time",1600,0,1600,1024,-256,256);
  s3List->Add(s3_rings_sectors_singlesTvT);
  s3_rings_sectors_singlesTvSec = new TH2F("s3_rings_sectors_singles_timing_v_sector","S3_sectorT_vs_ringT_Multiplicity_1_vs_Sector;Sector;Sector Time - Ring Time",32,0,32,1024,-256,256);
  s3List->Add(s3_rings_sectors_singlesTvSec);
  s3_rings_sectors_singlesTvRing = new TH2F("s3_rings_sectors_singles_timing_v_ring","S3_sectorT_vs_ringT_Multiplicity_1_vs_Ring;Ring;Sector Time - Ring Time",24,0,24,1024,-256,256);
  s3List->Add(s3_rings_sectors_singlesTvRing);
  s3_rings_sectors_singlesTvE = new TH2F("s3_rings_sectors_singles_timing_v_energy","S3_sectorT_vs_ringT_Multiplicity_1_vs_Energy;Energy (keV);Sector Time - Ring Time",4096,4096,8192,1024,-256,256);
  s3List->Add(s3_rings_sectors_singlesTvE);
  s3_rings_sectors_gated = new TH2F("s3_rings_sectors_gated","S3_sectorE_vs_ringE_gated_40keV;Ring Energy;Sector Energy",2048,0,8192,2048,0,8192);
  s3List->Add(s3_rings_sectors_gated);

  //TIGRESS-TIGRESS
  addT_addT = new TH1F("addT_addT","Tigress-Tigress_time",4096,-2048,2048); 
  tigtigList->Add(addT_addT);
  addE_addE = new TH2F("addE_addE", "Addback_Gamma_Gamma", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE);
  addE_addE_tg = new TH2F("addE_addE_tg", "Addback_Gamma_Gamma_Time_Gated", 4096, 0, 8192, 4096, 0, 8192);
  tigtigList->Add(addE_addE_tg);

  //TIGRESS-S3
  addT_s3T = new TH1F("addT_s3T","Tigress-S3_time",4096,-2048,2048); 
  tigs3List->Add(addT_s3T);
  excE = new TH1F("excE","Excitation_Energy_MeV",2800,-2.,12.); 
  tigs3List->Add(excE);
  excE_theta = new TH2F("excE_theta","Excitation_Energy_Theta",180,0,180,2800,-2.,12.); 
  tigs3List->Add(excE_theta);
  addDopp = new TH1F("addDopp", "Addback_Doppler_Energy", 8192, 0, 8192);
  tigs3List->Add(addDopp);
  addDopp_ANum = new TH2F("addDopp_ANum", "Addback_Doppler_Energy_V_ArrayNumber", 64, 0, 64, 8192, 0, 8192);
  tigs3List->Add(addDopp_ANum);
  addDopp_exc = new TH2F("addDopp_Exc", "Addback_Doppler_Energy_V_Excitation_Energy", 8192, 0, 8192, 2800, -2., 12.);
  tigs3List->Add(addDopp_exc);
  addE_s3_E  = new TH2F("addE_s3E", "Addback_Energy_V_S3_Energy", 8192, 0, 8192, 2048, 0, 32768);
  tigs3List->Add(addE_s3_E);
  addDopp_s3_E  = new TH2F("addDopp_s3E", "Addback_Doppler_Energy_V_S3_Energy", 8192, 0, 8192, 2048, 0, 32768);
  tigs3List->Add(addDopp_s3_E);
 

  // Raw EMMA
  xPos = new TH1F("xPos", "PGAC_X_Position", 160, -80, 80);
  emmaList->Add(xPos);
  yPos = new TH1F("yPos", "PGAC_Y_Position", 60, -30, 30);
  emmaList->Add(yPos);
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
    sprintf(hname, "EMMA_IC-Segment_%i", i);
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
    ssbE[i] = new TH1F(hname, hname, 16384, 0, 98304);
    ssbET[i] = new TH2F(hname2, hname2, 4000, 0, 4000, 16384, 0, 98304);
    emmaList->Add(ssbE[i]);
    emmaList->Add(ssbET[i]);
  }

  //ToF Spectra
  addemmatof = new TH1F("addemmatof", "Addback_Time - EMMA Time", 20000, -10000, 10000);
  tigemmaList->Add(addemmatof);
  s3emmatof = new TH1F("s3emmatof", "S3_Time - EMMA Time", 20000, -10000, 10000);
  tigemmaList->Add(s3emmatof);
  addE_tof = new TH2F("addE_tof", "Addback_Energy_vs_(Addback_Time - EMMA Time)", 20000, -10000, 10000, 8192, 0, 8192);
  tigemmaList->Add(addE_tof);
  tDE = new TH2F("tDE","tDE",4000,0,4000,400,-2000,2000);
  tigemmaList->Add(tDE);

  //Time Gated Spectra
  excE_tg = new TH1F("excE_tg","Excitation_Energy_MeV_Time_Gated",2800,-2.,12.); 
  tgList->Add(excE_tg);
  excE_theta_tg = new TH2F("excE_theta_tg","Excitation_Energy_Theta_Time_Gated",180,0,180,2800,-2.,12.); 
  tgList->Add(excE_theta_tg);
  addDopp_tg = new TH1F("addDopp_tg", "Addback_Doppler_Energy_Time_Gated", 8192, 0, 8192);
  tgList->Add(addDopp_tg);
  addDopp_ANum_tg = new TH2F("addDopp_ANum_tg", "Addback_Doppler_Energy_V_ArrayNumber_Time_Gated", 64, 0, 64, 8192, 0, 8192);
  tgList->Add(addDopp_ANum_tg);
  addDopp_exc_tg = new TH2F("addDopp_Exc_tg", "Addback_Doppler_Energy_V_Excitation_Energy_Time_Gated", 8192, 0, 8192, 2800, -2., 12.);
  tgList->Add(addDopp_exc_tg);
  addE_addE_tofg = new TH2F("addE_addE_tofg", "Addback_Gamma_Gamma_Time_Gated", 4096, 0, 8192, 4096, 0, 8192);
  tgList->Add(addE_addE_tofg);
  addDopp_addDopp_tg = new TH2F("addDopp_addDopp_tg", "Addback_Gamma_Gamma_Time_Gated_Doppler_Corrected", 4096, 0, 8192, 4096, 0, 8192);
  tgList->Add(addDopp_addDopp_tg);
  s3_E_theta_gated = new TH2F("s3_E_theta_gated","s3_E_theta_gated",180,0,180, 820, 0, 16384); 
  tgList->Add(s3_E_theta_gated);
  s3_E_gated = new TH1F("S3_E_gated","S3_E_gated",8192,0,16384);
  s3List->Add(s3_E); 
  hitmap_time_gated = new TH2F("hitmap_time_gated","Hitmap_gated_on_emma_tof",200,-50,50,200,-50,50);
  s3List->Add(hitmap_time_gated);
  icSumVSi_gated_350 = new TH2F("icSumVSi_gated_350","Si_Energy_VS_IC_SUM_gated_on_350_gamma",4096,0,4096,16384,0,16384);
  tgList->Add(icSumVSi_gated_350);



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

  //PID Gated Spectra
  excE_PIDG = new TH1F("excE_PIDG","Excitation_Energy_MeV_PID_Gated",2800,-2.,12.); 
  PIDgatedList->Add(excE_PIDG);
  excE_theta_PIDG = new TH2F("excE_theta_PIDG","Excitation_Energy_Theta_PID_Gated",180,0,180,2800,-2.,12.); 
  PIDgatedList->Add(excE_theta_PIDG);
  addDopp_PIDG = new TH1F("addDopp_PIDG", "Addback_Doppler_Energy_PID_Gated", 8192, 0, 8192);
  PIDgatedList->Add(addDopp_PIDG);
  addDopp_ANum_PIDG = new TH2F("addDopp_ANum_PIDG", "Addback_Doppler_Energy_V_ArrayNumber_PID_Gated", 64, 0, 64, 8192, 0, 8192);
  PIDgatedList->Add(addDopp_ANum_PIDG);
  addDopp_exc_PIDG = new TH2F("addDopp_exc_PIDG", "Addback_Doppler_Energy_V_Excitation_Energy_PID_Gated", 8192, 0, 8192, 2800, -2., 12.);
  PIDgatedList->Add(addDopp_exc_PIDG);

  excE_PIDG_bg = new TH1F("excE_PIDG_bg","Excitation_Energy_MeV_PID_Gated_Random_Time_Gated",2800,-2.,12.); 
  PIDgatedList->Add(excE_PIDG_bg);
  excE_theta_PIDG_bg = new TH2F("excE_theta_PIDG_bg","Excitation_Energy_Theta_PID_Gated_Random_Time_Gated",180,0,180,2800,-2.,12.); 
  PIDgatedList->Add(excE_theta_PIDG_bg);
  addDopp_PIDG_bg = new TH1F("addDopp_PIDG_bg", "Addback_Doppler_Energy_PID_Gated_Random_Time_Gated", 8192, 0, 8192);
  PIDgatedList->Add(addDopp_PIDG_bg);
  addDopp_ANum_PIDG_bg = new TH2F("addDopp_ANum_PIDG_bg", "Addback_Doppler_Energy_V_ArrayNumber_PID_Gated_Random_Time_Gated", 64, 0, 64, 8192, 0, 8192);
  PIDgatedList->Add(addDopp_ANum_PIDG_bg);
  addDopp_exc_PIDG_bg = new TH2F("addDopp_exc_PIDG_bg", "Addback_Doppler_Energy_V_Excitation_Energy_PID_Gated_Random_Time_Gated", 8192, 0, 8192, 2800, -2., 12.);
  PIDgatedList->Add(addDopp_exc_PIDG_bg);
  aE_icSumVSi_gated = new TH1F("aE_icSumVSi_gated","Addback_gamma_energy_gated_on_IC_Si_cut",8192,0,8192);
  PIDgatedList->Add(aE_icSumVSi_gated); 

  rings_PIDG = new TH2F("rings_PID_Gated","rings_PID_Gated",24,0,24,2048,0,8192); 
  PIDgatedListS3->Add(rings_PIDG);
  s3_E_theta_gated_PIDG = new TH2F("s3_E_theta_gated_PID_Gated","s3_E_theta_gated_PID_Gated",180,0,180, 820, 0, 16384); 
  PIDgatedListS3->Add(s3_E_theta_gated_PIDG);

}
