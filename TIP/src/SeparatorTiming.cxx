//Sort code to check TIP and TIGRESS timing windows
//windows are defined in common.h

#define SeparatorTiming_cxx
#include "common.h"
#include "SeparatorTiming.h"

using namespace std;

void SeparatorTiming::SortData(char const *afile, char const *calfile, char const *outfile){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  //setup the output file
  TFile *myfile = new TFile(outfile, "RECREATE");
  //setup the output tree
  TTree *outTree = AnalysisTree->GetTree()->CloneTree(0);
  //TTigress *outTig;
  //TTip *outTip 
  
  

  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }
  else
  {
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  TTip *tip = 0;
  if (AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
  }

  //TBranch *outTigBranch = outTree->Branch("TTigress",&tigress);
  //TBranch *outTipBranch = outTree->Branch("TTip",&tip);

  bool suppAdd = false;

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  printf("\nSorting analysis events...\n");
  for (int jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    if(tigress && tip){

      if(tip->GetMultiplicity()>MAXNUMTIPHIT){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }
      if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }

      bool pass = false;
      uint64_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup
      if(passedtimeGate&(1ULL<<61)){
        if(passedtimeGate&(1ULL<<62)){
          if(passedtimeGate&(1ULL<<63)){
            pass = true;
          }
        }
      }

      if(pass){

        /*for(int tigHitIndAB = 0; tigHitIndAB < MAXNUMTIGHIT; tigHitIndAB++){
          if(tigHitIndAB < tigress->GetAddbackMultiplicity()){
            if(!(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT)))){
              //cout << "num tig hits before: " << tigress->GetAddbackMultiplicity() << endl;
              tigress->GetAddbackHit(tigHitIndAB)->Clear();
              //cout << "num tig hits after: " << tigress->GetAddbackMultiplicity() << endl;
            }
          }
        }*/

        for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
          if(!(passedtimeGate&(1ULL<<tipHitInd))){
            tip->GetTipHit(tipHitInd)->Clear();
            //tip->GetHitVector().clear(tip->GetHitVector().begin()+tipHitInd);
          }
        }

        tigress->ResetAddback(); //needed so that addback is reconstructed in subsequent analysis

        outTree->Fill();
      }

    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  cout << "Entries in separated output data: " << outTree->GetEntries() << endl;

  cout << "Writing histograms to " << outfile << endl;

  myfile->Write();
  myfile->Close();

}
int main(int argc, char **argv){

  SeparatorTiming *mysort = new SeparatorTiming();

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0)
  {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1)
  {
    cout << "Arguments: SeparatorTiming analysis_tree calibration_file output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }
  else if (argc == 2)
  {
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }
  else
  {
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
