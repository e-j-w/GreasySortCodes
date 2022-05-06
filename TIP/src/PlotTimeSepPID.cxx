//Sort code to pot TIP PID for time separated data
//timing windows are defined in common.h

#define PlotTimeSepPID_cxx
#include "common.h"
#include "PlotTimeSepPID.h"

using namespace std;

float lastTIPHitT[NTIP]; //stores the last hit time for each detector
Int_t numTipRingPileupHits[NTIPRING], numTipRingHits[NTIPRING];


void PlotTimeSepPID::SortData(char const *afile, char const *calfile, char const *outfile){

  Initialise();

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen())
  {
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

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

  unsigned long int numTipHits = 0;
  unsigned long int numTipPileupHits = 0;

  double_t tipPID = -1000.0;

  //Defining Pointers
  TTipHit *tip_hit;
  const std::vector<Short_t> *wf; //for CsI waveform

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

      uint32_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup

      for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
        if(passedtimeGate&(1U<<tipHitInd)){
          tip_hit = tip->GetTipHit(tipHitInd);
          numTipHits++;

          wf = tip_hit->GetWaveform();
          TPulseAnalyzer pulse;
          pulse.SetData(*wf, 0);
          if(wf->size() > 50){
            tipPID = pulse.CsIPID();
            if(tipPID > -1000.0) //in (modified) GRSISort, failed fits given value of -1000.0
              tipPID += 100.;
            
          }

          //PID related stuff
          if(tipPID>=0){ //PID was found
            tip_E_PID_Sum->Fill(tip_hit->GetEnergy(),tipPID);
            if((tip_hit->GetTipChannel() > 0)&&(tip_hit->GetTipChannel() <= NTIP)){
              tip_E_PID_Ring[getTIPRing(tip_hit->GetTipChannel())]->Fill(tip_hit->GetEnergy(),tipPID);
              tip_E_PID[tip_hit->GetTipChannel()-1]->Fill(tip_hit->GetEnergy(),tipPID);
            }
          }
        }
        
      }


    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP pileup hits: " << numTipPileupHits << " (" << (float)(numTipPileupHits*100)/float(numTipPileupHits + numTipHits) << " %)" << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();

  TDirectory *tippiddir = myfile->mkdir("TIP PID");
  tippiddir->cd();
  tipPIDList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
}
int main(int argc, char **argv){

  PlotTimeSepPID *mysort = new PlotTimeSepPID();

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
    cout << "Arguments: PlotTimeSepPID analysis_tree calibration_file output_file" << endl;
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
