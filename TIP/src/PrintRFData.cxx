//Sort code to draw RF waveforms event-by-event

#include <iostream>
#include <iomanip>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TRF.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

#define PrintRF_cxx

using namespace std;

TApplication *theApp;

void SortData(char const *afile, char const *calfile){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if(!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long analentries = AnalysisTree->GetEntries();

  TRF *rf = 0;
  if(AnalysisTree->FindBranch("TRF")){
    AnalysisTree->SetBranchAddress("TRF", &rf);
  }else{
    cout << "Branch 'TRF' not found! TRF variable is NULL pointer" << endl;
  }

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  printf("\nSorting analysis events...\n");
  for(Long64_t jentry = 0; jentry < analentries; jentry++){

    AnalysisTree->GetEntry(jentry);
    
    if(rf){
      cout << "Entry " << jentry << endl;
      cout << "ts: " << rf->TimeStamp() << ", RF period: " << rf->Period() << ", Fit time: " << rf->Time() << " ns." << endl;
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "No more events in this run." << endl;

}
int main(int argc, char **argv)
{

  char const *afile;
  char const *calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if(argc == 1){
    cout << "Arguments: PrintRFData analysis_tree calibration_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else if(argc == 2){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    printf("Analysis file: %s\nCalibration file: %s\n", afile, calfile);
  }else if(argc == 3){
    afile = argv[1];
    calfile = argv[2];
    printf("Analysis file: %s\nCalibration file: %s\n", afile, calfile);
  }else{
    printf("Arguments: PrintRFData analysis_tree calibration_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  SortData(afile, calfile);

  return 0;
}
