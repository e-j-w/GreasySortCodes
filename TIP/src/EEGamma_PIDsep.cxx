//Plots TIGRESS gamma-gamma matrices for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_PIDsep_cxx
#include "common.h"
#include "EEGamma_PIDsep.h"

using namespace std;

void EEGamma_PIDsep::SortData(char const *afile, char const *calfile, char const *outfile){

  Initialise();

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if(!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  TTip *tip = 0;
  if(AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
  }

  unsigned long int numTipHits = 0;
  unsigned long int numTigABHits = 0;

  unsigned int evtNumProtons, evtNumAlphas;

  //Defining Pointers
  TTipHit *tip_hit;
  TTigressHit *add_hit, *add_hit2;

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  printf("\nSorting analysis events...\n");
  for(Long64_t jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    evtNumProtons = 0;
    evtNumAlphas = 0;

    if(tigress && tip){

      if(tip->GetMultiplicity()>MAXNUMTIPHIT){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }
      if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }

      uint64_t passedtimeGate = passesTimeGate(tigress,tip,2,2); //also rejects pileup

      //count the number of protons or alphas
      for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
        if(passedtimeGate&(1ULL<<tipHitInd)){
          tip_hit = tip->GetTipHit(tipHitInd);
          numTipHits++;

          switch(getParticleType(tip_hit,gates)){ //see common.cxx
            case 4:
              evtNumAlphas++;
              break;
            case 1:
              evtNumProtons++;
              break;
            case 0:
            default:
              break;
          }

        }
      }
      
      if(evtNumProtons<=MAX_NUM_PARTICLE){
        if(evtNumAlphas<=MAX_NUM_PARTICLE){
          if(evtNumProtons+evtNumAlphas<=MAX_NUM_PARTICLE){

            for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                numTigABHits++;
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                if(!add_hit->BGOFired() && add_hit->GetEnergy() > MIN_TIG_EAB){
                  //TIGRESS PID separated addback energy
                  for(int tigHitIndAB2=tigHitIndAB+1;tigHitIndAB2<tigress->GetAddbackMultiplicity();tigHitIndAB2++){
                    if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                      add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                      if(!add_hit2->BGOFired() && add_hit2->GetEnergy() > MIN_TIG_EAB){
                        tigEE_xayp[evtNumProtons][evtNumAlphas]->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                        tigEE_xayp[evtNumProtons][evtNumAlphas]->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
                      }
                    }
                  }
                }
              }
            }
          }else{
            cout << "Event " << jentry << " has too many charged particles (" << evtNumProtons+evtNumAlphas << ")!" << endl;
          }
        }else{
          cout << "Event " << jentry << " has too many alphas (" << evtNumAlphas << ")!" << endl;
        }
      }else{
        cout << "Event " << jentry << " has too many protons (" << evtNumProtons << ")!" << endl;
      }
        
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIGRESS addback hits: " << numTigABHits << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();

  TDirectory *tigpidsepdir = myfile->mkdir("TIGRESS PID Separated");
  tigpidsepdir->cd();
  tigPIDSepList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
  analysisfile->Close();
}

int main(int argc, char **argv){

  EEGamma_PIDsep *mysort = new EEGamma_PIDsep();

  char const *afile;
  char const *outfile;
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
    cout << "Plots TIGRESS spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_PIDsep analysis_tree calibration_file output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else if(argc == 2){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if(argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if(argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else{
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
