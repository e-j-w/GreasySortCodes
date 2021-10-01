//g++ clock_check.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o Clock

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TEmma.h"
#include "TGRSIDetectorHit.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

void clock_check(char const * infile, char const * calfile, char const * outfile) {

 
  TList * list = new TList;

  TH1F *tD = new TH1F("tD","tD",4096,-2048,2048);list->Add(tD);
  TH2F *tDE = new TH2F("tDE","tDE",4000,0,4000,1000,-5000,5000);list->Add(tDE);
  TH1F *tTig = new TH1F("EMMA trigger time at TIGRESS","EMMA trigger time at TIGRESS",4096,0,2048);
  tTig->GetYaxis()->SetTitle("Counts");tTig->GetXaxis()->SetTitle("Time (s)");list->Add(tTig);
  TH1F *tAnode = new TH1F("EMMA anode hit time","EMMA anode hit time",4096,0,2048);
  tAnode->GetYaxis()->SetTitle("Counts");tAnode->GetXaxis()->SetTitle("Time (s)");list->Add(tAnode);
  TH1F *tIC = new TH1F("EMMA IC hit time","EMMA IC hit time",4096,0,2048);
  tIC->GetYaxis()->SetTitle("Counts");tAnode->GetXaxis()->SetTitle("Time (s)");list->Add(tIC);
  TH1F *tSi = new TH1F("EMMA Si hit time","EMMA Si hit time",4096,0,2048);
  tSi->GetYaxis()->SetTitle("Counts");tAnode->GetXaxis()->SetTitle("Time (s)");list->Add(tSi);

  tD->GetXaxis()->SetTitle("TIGRESS DAQ Time -- EMMA Daq Time ns");
  tDE->GetXaxis()->SetTitle("Run Time");
  tDE->GetYaxis()->SetTitle("TIGRESS DAQ Time -- EMMA Daq Time ns");
  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    printf("Opening file failed, aborting\n");
    return;
  }
  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  printf("%i tree files, details:\n", AnalysisTree->GetNtrees());
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  uint64_t numAnodeHits = 0;
  uint64_t numSiHits = 0;
  uint64_t numICHits = 0;
  uint64_t numTriggerHits = 0;

  TEmma * emma = nullptr;
  if (AnalysisTree->FindBranch("TEmma")) {
    AnalysisTree->SetBranchAddress("TEmma", & emma);
  } else {
    cout << "Branch 'TEmma' not found! TEmma variable is NULL pointer!" << endl;
  }

  int i=0;
  TEmmaHit * trigger_hit;
  TEmmaHit * anode_hit, *ichit, *si_hit, *ic_hit;
  printf("Begin sort\n");
  for (int jentry = 0; jentry < (nentries - 1); jentry++) {
    tree->GetEntry(jentry);
    for (int one = 0; one < emma->GetTriggerMultiplicity(); one++) {
      trigger_hit = emma->GetTriggerHit(one);
      //cout << "Trigger hit " << one << ", energy: " << trigger_hit->GetEnergy() << endl;
      if(trigger_hit->GetEnergy() > 2000)continue;
      numTriggerHits++;
      tTig->Fill(trigger_hit->GetTime()/pow(10,9)); //cout << "filling " << trigger_hit->GetTime()/pow(10,9) << endl;
      for (int two = 0; two < emma->GetTdcMultiplicity(); two++) {
        anode_hit = emma->GetTdcHit(two);
        //cout << "Anode hit segment: " << anode_hit->GetSegment() << endl;
        //if(anode_hit->GetSegment()==2){
          //cout << "Anode hit " << two << ", segment: " << anode_hit->GetSegment() << endl;
          //cout << "times: " << anode_hit->GetTime() << ", " << trigger_hit->GetTime() << ", diff: " << anode_hit->GetTimeStamp()*anode_hit->GetTimeStampUnit()-trigger_hit->GetTimeStamp()*trigger_hit->GetTimeStampUnit() << endl;
          tD->Fill(anode_hit->GetTime()-trigger_hit->GetTime());
          tDE->Fill(trigger_hit->GetTime()/pow(10,9),anode_hit->GetTimeStamp()*anode_hit->GetTimeStampUnit()-trigger_hit->GetTimeStamp()*trigger_hit->GetTimeStampUnit());
        //}
      }
    }
    for (int two = 0; two < emma->GetTdcMultiplicity(); two++) {
      anode_hit = emma->GetTdcHit(two);
      //cout << "Anode hit segment: " << anode_hit->GetSegment() << endl;
      //if(anode_hit->GetSegment()==2){
        tAnode->Fill(anode_hit->GetTimeStamp()*anode_hit->GetTimeStampUnit()/pow(10,9)); //cout << "filling anode " << anode_hit->GetTimeStamp()*anode_hit->GetTimeStampUnit()/pow(10,9) << endl;
      //}
      /*for (int k = 0; k< emma->GetICMultiplicity(); k++) {
        anode_hit = emma->GetTdcHit(two);
        ichit = emma->GetICHit(two);
        //cout << "tdc hit " << two << " time: " << anode_hit->GetTime() << ", ic hit " << k << " time: " << ichit->GetTime() << ", diff: " << anode_hit->GetTime()-ichit->GetTime() << endl;
      }*/
      numAnodeHits++;
    }
    for(int e=0;e<emma->GetSiMultiplicity();e++){
	    si_hit = emma->GetSiHit(e);
      tSi->Fill(si_hit->GetTimeStamp()*si_hit->GetTimeStampUnit()/pow(10,9));
      numSiHits++;
    }
    for (int j = 0; j < emma->GetICMultiplicity(); j++){
      ic_hit = emma->GetICHit(j);
      tIC->Fill(ic_hit->GetTimeStamp()*ic_hit->GetTimeStampUnit()/pow(10,9));
      numICHits++;
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "TIGRESS Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }
  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";
  cout << "Histograms written, sorting complete" << endl;

  cout << "Number of trigger hits: " << numTriggerHits << endl;
  cout << "Number of anode hits: " << numAnodeHits << endl;
  cout << "Number of IC hits: " << numICHits << endl;
  cout << "Number of focal plane Si hits: " << numSiHits << endl;

  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();
}

int main(int argc, char ** argv) {

  char const * afile;
  char const * outfile;
  char const * calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0) {
	  grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
	  cout << "Insufficient arguments, provide argument tree files" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
	  calfile = "CalibrationFile.cal";
	  outfile = "Clock.root";
  } else if (argc == 3) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = "Clock.root";
  } else if (argc == 4) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = argv[3];
  } else if (argc > 4) {
	  printf("Too many arguments\n");
	  return 0;
  }

  printf("Input file:%s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

  TParserLibrary::Get()->Load();

  clock_check(afile, calfile, outfile);

  return 0;
}
