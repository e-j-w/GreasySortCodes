//Plots TIGRESS gamma-gamma matrices for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define TTGamma_PIDsep_cxx
#include "common.h"
#include "TTGamma_PIDsep.h"

using namespace std;

void TTGamma_PIDsep::Initialise(const int nP, const int nA, const float Elow_start, const float Ehigh_start, const float Elow_stop, const float Ehigh_stop){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tttimeList = new TList;

  tigTT = new TH1F(Form("TIGRESS-TIGRESS timing (%ip%ia gate, start %8.2f-%8.2f keV, stop  %8.2f-%8.2f keV)",nP,nA,Elow_start,Ehigh_start,Elow_stop,Ehigh_stop),Form("TIGRESS-TIGRESS timing (%ip%ia gate, start %8.2f-%8.2f keV, stop  %8.2f-%8.2f keV)",nP,nA,Elow_start,Ehigh_start,Elow_stop,Ehigh_stop),4096,-1024,1024);
  tigTT->GetXaxis()->SetTitle("t_{diff} (ns)");
  tigTT->GetYaxis()->SetTitle("Counts");
  tttimeList->Add(tigTT);

  //Setup TIP PID gates
  gates = new PIDGates;
}

void TTGamma_PIDsep::WriteData(char const *outfile){

  cout << "Event sorting complete" << endl;
  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();

  TDirectory *tttimedir = myfile->mkdir("TIG-TIG timing");
  tttimedir->cd();
  tttimeList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();

}

void TTGamma_PIDsep::SortData(char const *afile, char const *calfile, const int nP, const int nA, const float Elow_start, const float Ehigh_start, const float Elow_stop, const float Ehigh_stop){

  printf("Number of protons: %i, alphas: %i\n", nP, nA);
  printf("Start gate: %8.2f - %8.2f keV\n", Elow_start, Ehigh_start);
  printf("Stop gate:  %8.2f - %8.2f keV\n", Elow_stop, Ehigh_stop);

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if(!analysisfile->IsOpen()){
    cout << "Opening file " << afile << " failed, aborting" << endl;
    return;
  }

  cout << "File " << afile << " opened" << endl;
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if(AnalysisTree->FindBranch("TTigress")){
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

  unsigned int evtNumProtons, evtNumAlphas;

  //Defining Pointers
  TTipHit *tip_hit;
  TTigressHit *add_hit, *add_hit2;

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  Long64_t numSeparatedEvents = 0;

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

      uint64_t passedtimeGate = passesTimeGateAB(tigress,tip,2,2); //also rejects pileup

      if(passedtimeGate&(1ULL<<TIPTIGFLAG)){
        numSeparatedEvents++;
        //count the number of protons or alphas
        for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
          if(passedtimeGate&(1ULL<<tipHitInd)){
            tip_hit = tip->GetTipHit(tipHitInd);

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
        
        if(evtNumProtons == nP){
          if(evtNumAlphas == nA){
            if(evtNumProtons<=MAX_NUM_PARTICLE){
              if(evtNumAlphas<=MAX_NUM_PARTICLE){
                if(evtNumProtons+evtNumAlphas<=MAX_NUM_PARTICLE){
                  for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
                    if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                      add_hit = tigress->GetAddbackHit(tigHitIndAB);
                      //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                      if(!add_hit->BGOFired() && add_hit->GetEnergy() > MIN_TIG_EAB){
                        if(add_hit->GetEnergy() >= Elow_start){
                          if(add_hit->GetEnergy() <= Ehigh_start){
                            for(int tigHitIndAB2=0;tigHitIndAB2<tigress->GetAddbackMultiplicity();tigHitIndAB2++){
                              if(tigHitIndAB2!=tigHitIndAB){
                                if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                                  add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                                  if(!add_hit2->BGOFired() && add_hit2->GetEnergy() > MIN_TIG_EAB){
                                    if(add_hit2->GetEnergy() >= Elow_stop){
                                      if(add_hit2->GetEnergy() <= Ehigh_stop){
                                        cout << "Hit 1: " << add_hit->GetEnergy() << " keV, t=" << add_hit->GetTime() << ", Hit 2: " << add_hit2->GetEnergy() << " keV, t=" << add_hit2->GetTime() << ", tdiff=" << (add_hit2->GetTime() - add_hit->GetTime()) << endl;
                                        tigTT->Fill(add_hit2->GetTime() - add_hit->GetTime());
                                      }
                                    }
                                  }
                                }
                              }
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
        }
      }
      
        
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of time separated events: " << numSeparatedEvents << endl;

  analysisfile->Close();
  
}

int main(int argc, char **argv){

  TTGamma_PIDsep *mysort = new TTGamma_PIDsep();

  char const *afile;
  char const *outfile;
  char const *calfile;
  float Elow_start = 0.0f;
  float Ehigh_start = 0.0f;
  float Elow_stop = 0.0f;
  float Ehigh_stop = 0.0f;
  int nP = 0;
  int nA = 0;
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
    cout << "Plots TIGRESS timing spectra for PID and time separated data." << endl;
    cout << "Arguments: TTGamma_PIDsep analysis_tree calibration_file nP nA Elow_start EHigh_start Elow_stop Ehigh_stop output_file" << endl;
    cout << "  *analysis_tree* can be a single tree (extension .root) or a" << endl;
    cout << "  plaintext list of trees (one per line, extension .list)." << endl;
    return 0;
  }else if(argc == 9){
    afile = argv[1];
    calfile = argv[2];
    nP = atoi(argv[3]);
    nA = atoi(argv[4]);
    Elow_start = (float)atof(argv[5]);
    Ehigh_start = (float)atof(argv[6]);
    Elow_stop = (float)atof(argv[7]);
    Ehigh_stop = (float)atof(argv[8]);
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if(argc == 10){
    afile = argv[1];
    calfile = argv[2];
    nP = atoi(argv[3]);
    nA = atoi(argv[4]);
    Elow_start = (float)atof(argv[5]);
    Ehigh_start = (float)atof(argv[6]);
    Elow_stop = (float)atof(argv[7]);
    Ehigh_stop = (float)atof(argv[8]);
    outfile = argv[9];
  }else{
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  if((nP < 0)||(nA < 0)){
    printf("ERROR: invalid number of protons or alphas.\n");
    return 0;
  }
  if(Elow_start > Ehigh_start){
    printf("ERROR: start gate is invalid.\n");
    return 0;
  }
  if(Elow_stop > Ehigh_stop){
    printf("ERROR: stop gate is invalid.\n");
    return 0;
  }

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(strcmp(dot + 1, "root") == 0){

    //single analysis tree
    cout << "Analysis tree: " << afile << endl;
    printf("Calibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

    theApp=new TApplication("App", &argc, argv);
    mysort->Initialise(nP, nA, Elow_start, Ehigh_start, Elow_stop, Ehigh_stop);

    mysort->SortData(afile, calfile, nP, nA, Elow_start, Ehigh_start, Elow_stop, Ehigh_stop);

  }else if(strcmp(dot + 1, "list") == 0){

    //list of analysis trees
    cout << "Analysis tree list: " << afile << endl;
    printf("Calibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

    theApp=new TApplication("App", &argc, argv);
    mysort->Initialise(nP, nA, Elow_start, Ehigh_start, Elow_stop, Ehigh_stop);

    FILE *listfile;
    char str[256];

    if((listfile=fopen(afile,"r"))==NULL){
      cout << "ERROR: Cannot open the input file: " << afile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          mysort->SortData(str, calfile, nP, nA, Elow_start, Ehigh_start, Elow_stop, Ehigh_stop);
        }
      }
    }

  }else{
    cout << "ERROR: Invalid analysis tree file type!" << endl;
    return 0;
  }

  mysort->WriteData(outfile);

  return 0;
}
