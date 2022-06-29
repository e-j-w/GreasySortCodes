//Sort code to check TIP and TIGRESS timing windows
//windows are defined in common.h

#define CheckTimingWindows_cxx
#include "common.h"
#include "CheckTimingWindows.h"

using namespace std;

void CheckTimingWindows::SortData(char const *afile, char const *calfile, char const *outfile){

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
  unsigned long int numTigABHits = 0;
  unsigned long int numTipTigHits = 0;
  unsigned long int numPassTipTip = 0;
  unsigned long int numPassTigTig = 0;
  unsigned long int numPassTipTig = 0;
  Double_t tipFitTimes[MAXNUMTIPHIT];

  //Defining Pointers
  TTigressHit *add_hit, *add_hit2;
  TTipHit *tip_hit, *tip_hit2;
  
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

      numTipTigHits++;

      uint64_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup
      if(passedtimeGate&(1ULL<<61)){
        numPassTipTip++;
      }
      if(passedtimeGate&(1ULL<<62)){
        numPassTigTig++;
      }
      if(passedtimeGate&(1ULL<<63)){
        numPassTipTig++;
      }

      //count the number of hits which passed the timing gates
      Int_t tipMultPassed=0;
      Int_t tigMultPassed=0;
      Int_t tigSuppMult=0;
      Int_t tigSuppMultPassed=0;
      for(int tipHitInd = 0; tipHitInd < MAXNUMTIPHIT; tipHitInd++){
        if(tipHitInd < tip->GetMultiplicity()){
          if(passedtimeGate&(1ULL<<tipHitInd)){
            tipMultPassed++;
          }
        }
      }
      for(int tigHitIndAB = 0; tigHitIndAB < MAXNUMTIGHIT; tigHitIndAB++){
        if(tigHitIndAB < tigress->GetAddbackMultiplicity()){
          if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
            tigMultPassed++;
          }
          add_hit = tigress->GetAddbackHit(tigHitIndAB);
          suppAdd = add_hit->BGOFired();
          //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
          if (!suppAdd && add_hit->GetEnergy() > 15){
            tigSuppMult++;
            if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
              tigSuppMultPassed++;
            }
          }
        }
      }

      //get fit times
      for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
        tip_hit = tip->GetTipHit(tipHitInd);
        tipFitTimes[tipHitInd] = getTipFitTime(tip_hit,tip_waveform_pretrigger);
      }

      for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
        numTipHits++;
        tip_hit = tip->GetTipHit(tipHitInd);
        for(int tipHitInd2 = tipHitInd+1; tipHitInd2 < tip->GetMultiplicity(); tipHitInd2++){
          tip_hit2 = tip->GetTipHit(tipHitInd2);
          Double_t tDiff = tipFitTimes[tipHitInd] - tipFitTimes[tipHitInd2];
          tiptipFitT->Fill(tDiff);
          //tiptipFitT->Fill(tip_hit->GetTime() - tip_hit2->GetTime());
          if(passedtimeGate&(1ULL<<tipHitInd)){
            if(passedtimeGate&(1ULL<<tipHitInd2)){
              tiptipFitTPassed->Fill(tDiff);
            }
          }
        }
      }
      
      for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
        numTigABHits++;
        add_hit = tigress->GetAddbackHit(tigHitIndAB);
        suppAdd = add_hit->BGOFired();
        //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
        if (!suppAdd && add_hit->GetEnergy() > 15){

          //check time-correlated TIP events
          for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
            tip_hit = tip->GetTipHit(tipHitInd);
            Double_t tDiff = tipFitTimes[tipHitInd] - add_hit->GetTime();
            tipT_tigT_diff->Fill(tDiff);
            if(passedtimeGate&(1ULL<<tipHitInd)){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                tipT_tigT_diffPassed->Fill(tDiff);
              }
            }
            
          }

          //TIGRESS-TIGRESS addback
          for (int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < tigress->GetAddbackMultiplicity(); tigHitIndAB2++){

            add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
            suppAdd = add_hit2->BGOFired();
            if (!suppAdd && add_hit2->GetEnergy() > 15){
              Double_t tDiff = add_hit->GetTime() - add_hit2->GetTime();
              addT_addT->Fill(tDiff);
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                  addT_addTPassed->Fill(tDiff);
                }
              }
            }
          }
          
        }
      }

      tiptig_mult->Fill(tip->GetMultiplicity(),tigress->GetMultiplicity());
      tiptig_multSupp->Fill(tip->GetMultiplicity(),tigSuppMult);
      tiptig_multPassed->Fill(tipMultPassed,tigMultPassed);
      tiptig_multSuppPassed->Fill(tipMultPassed,tigSuppMultPassed);

    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigABHits << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP + TIGRESS events: " << numTipTigHits << endl << endl;
  cout << "Events passing TIP-TIP time gate: " << numPassTipTip << " (" << 100.0*numPassTipTip/(1.0*numTipTigHits) << " %)" << endl;
  cout << "Events passing TIG-TIG time gate: " << numPassTigTig << " (" << 100.0*numPassTigTig/(1.0*numTipTigHits) << " %)" << endl;
  cout << "Events passing TIP-TIG time gate: " << numPassTipTig << " (" << 100.0*numPassTipTig/(1.0*numTipTigHits) << " %)" << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();

  TDirectory *tiptipdir = myfile->mkdir("TIP-TIP");
  tiptipdir->cd();
  tiptipList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  TDirectory *tiptigdir = myfile->mkdir("TIP-TIGRESS");
  tiptigdir->cd();
  tiptigList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
}
int main(int argc, char **argv){

  CheckTimingWindows *mysort = new CheckTimingWindows();

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
    cout << "Arguments: CheckTimingWindows analysis_tree calibration_file output_file" << endl;
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
