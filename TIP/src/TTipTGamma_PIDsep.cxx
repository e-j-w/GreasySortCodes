//Plots TIGRESS gamma-gamma matrices for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define TTipTGamma_PIDsep_cxx
#include "common.h"
#include "TTipTGamma_PIDsep.h"

using namespace std;

void TTipTGamma_PIDsep::Initialise(const int nP, const int nA){

  printf("Start initialization\n");
  printf("Creating lists\n");

  tttimeList = new TList;

  tipTtigT_tigE = new TH2F(Form("TIP-TIGRESS timing (%ip%ia gate) vs. Addback energy",nP,nA),Form("TIP-TIGRESS timing (%ip%ia gate) vs. Addback energy",nP,nA),8192,0,8192,1024,896,1920);
  tipTtigT_tigE->GetYaxis()->SetTitle("t_{TIG} - t_{TIP} (ns)");
  tipTtigT_tigE->GetXaxis()->SetTitle("Addback energy (keV)");
  tttimeList->Add(tipTtigT_tigE);

  tipTtigT_doppE = new TH2F(Form("TIP-TIGRESS timing (%ip%ia gate) vs. Doppler energy",nP,nA),Form("TIP-TIGRESS timing (%ip%ia gate) vs. Doppler energy",nP,nA),8192,0,8192,1024,896,1920);
  tipTtigT_doppE->GetYaxis()->SetTitle("t_{TIG} - t_{TIP} (ns)");
  tipTtigT_doppE->GetXaxis()->SetTitle("Doppler corrected addback energy (keV)");
  tttimeList->Add(tipTtigT_doppE);

  EEGamma = new TH2F("Addback energy","Addback energy",4096,0,8192,4096,0,8192);
  EEGamma->GetYaxis()->SetTitle("Addback energy (keV)");
  EEGamma->GetXaxis()->SetTitle("Addback energy (keV)");
  tttimeList->Add(EEGamma);

  EEDopp = new TH2F("Doppler energy vs. Doppler energy","Doppler energy vs. Doppler energy",4096,0,8192,4096,0,8192);
  EEDopp->GetYaxis()->SetTitle("Doppler corrected addback energy (keV)");
  EEDopp->GetXaxis()->SetTitle("Doppler corrected addback energy (keV)");
  tttimeList->Add(EEDopp);

  //Setup TIP PID gates
  gates = new PIDGates;
}

void TTipTGamma_PIDsep::WriteData(char const *outfile){

  cout << "Event sorting complete" << endl;
  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();

  TDirectory *tttimedir = myfile->mkdir("TIP-TIG timing");
  tttimedir->cd();
  tttimeList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();

}

void TTipTGamma_PIDsep::SortData(char const *afile, char const *calfile, const int nP, const int nA){

  printf("Number of protons: %i, alphas: %i\n", nP, nA);

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    cout << "Opening file " << afile << " failed, aborting" << endl;
    return;
  }

  cout << "File " << afile << " opened" << endl;
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  TTip *tip = 0;
  if (AnalysisTree->FindBranch("TTip")){
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

      uint64_t passedtimeGate = passesTimeGate(tigress,tip,1,(uint8_t)(nP+nA)); //also rejects pileup

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
        
        double startTime = -1e30;
        if(evtNumProtons == nP){
        //if((evtNumProtons >= 1) && (evtNumProtons <= 2)){
          if(evtNumAlphas == nA){
            if(evtNumProtons<=MAX_NUM_PARTICLE){
              if(evtNumAlphas<=MAX_NUM_PARTICLE){
                if(evtNumProtons+evtNumAlphas<=MAX_NUM_PARTICLE){
                  for(int tipHitInd=0;tipHitInd<tigress->GetAddbackMultiplicity();tipHitInd++){
                    if(passedtimeGate&(1ULL<<tipHitInd)){
                      double fitTime = getTipFitTime(tip->GetTipHit(tipHitInd),tip_waveform_pretrigger);
                      if(fitTime > startTime){
                        startTime = fitTime; //the TIP time is the time of the last particle (closest in time to the gamma decay)
                      }
                    }
                  }
                  for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
                    if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                      add_hit = tigress->GetAddbackHit(tigHitIndAB);
                      //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                      if(!add_hit->BGOFired() && add_hit->GetEnergy() > 15){
                        tipTtigT_tigE->Fill(add_hit->GetEnergy(),add_hit->GetTime() - startTime);
                        tipTtigT_doppE->Fill(getEDoppFusEvap(add_hit,tip,passedtimeGate,gates),add_hit->GetTime() - startTime);
                      }
                    }
                  }
                  if(tigress->GetAddbackMultiplicity()>=2){
                    for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
                      if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                        add_hit = tigress->GetAddbackHit(tigHitIndAB);
                        //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                        if(!add_hit->BGOFired() && add_hit->GetEnergy() > 15){
                          for(int tigHitIndAB2=tigHitIndAB+1;tigHitIndAB2<tigress->GetAddbackMultiplicity();tigHitIndAB2++){
                            if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                              add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                              if(!add_hit2->BGOFired() && add_hit2->GetEnergy() > 15){
                                EEGamma->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                                EEDopp->Fill(getEDoppFusEvap(add_hit,tip,passedtimeGate,gates),getEDoppFusEvap(add_hit2,tip,passedtimeGate,gates));
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

  TTipTGamma_PIDsep *mysort = new TTipTGamma_PIDsep();

  char const *afile;
  char const *outfile;
  char const *calfile;
  int nP = 0;
  int nA = 0;
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
  if(argc == 1){
    cout << "Plots TIP-TIGRESS timing spectra for PID and time separated data." << endl;
    cout << "Arguments: TTipTGamma_PIDsep analysis_tree calibration_file nP nA output_file" << endl;
    cout << "  *analysis_tree* can be a single tree (extension .root) or a" << endl;
    cout << "  plaintext list of trees (one per line, extension .list)." << endl;
    return 0;
  }else if(argc == 5){
    afile = argv[1];
    calfile = argv[2];
    nP = atoi(argv[3]);
    nA = atoi(argv[4]);
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if(argc == 6){
    afile = argv[1];
    calfile = argv[2];
    nP = atoi(argv[3]);
    nA = atoi(argv[4]);
    outfile = argv[5];
  }else{
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  if((nP < 0)||(nA < 0)){
    cout << "ERROR: invalid number of protons or alphas." << endl;
    return 0;
  }

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(strcmp(dot + 1, "root") == 0){

    //single analysis tree
    cout << "Analysis tree: " << afile << endl;
    printf("Calibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

    theApp=new TApplication("App", &argc, argv);
    mysort->Initialise(nP, nA);

    mysort->SortData(afile, calfile, nP, nA);

  }else if(strcmp(dot + 1, "list") == 0){

    //list of analysis trees
    cout << "Analysis tree list: " << afile << endl;
    printf("Calibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

    theApp=new TApplication("App", &argc, argv);
    mysort->Initialise(nP, nA);

    FILE *listfile;
    char str[256];

    if((listfile=fopen(afile,"r"))==NULL){
      cout << "ERROR: Cannot open the input file: " << afile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          mysort->SortData(str, calfile, nP, nA);
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
