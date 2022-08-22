//Generates TIGRESS gamma-gated gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EGamma_PIDsep_mca_cxx
#include "common.h"
#include "EGamma_PIDsep_mca.h"

using namespace std;

void EGamma_PIDsep_mca::WriteData(const unsigned int numP, const unsigned int numA, const int gateEmin, const int gateEmax, const int writeProj){

  char const *outName = Form("EGamma_%up%ua_c%iw%i_gated.fmca",numP,numA,(gateEmax+gateEmin)/2,gateEmax-gateEmin);
  cout << "Writing gated histogram to: " << outName << endl;

  FILE *out;
  if((out = fopen(outName, "w")) == NULL){ //open the file
    cout << "ERROR: Cannot open the output file: " << outName << endl;
    return;
  }else{
    fwrite(&mcaOut,sizeof(mcaOut),1,out);
    fclose(out);
  }

  if(writeProj){

    char const *projName = Form("EGamma_%up%ua_proj.fmca",numP,numA);
    cout << "Writing projection histogram to: " << projName << endl;

    FILE *projOut;
    if((projOut = fopen(projName, "w")) == NULL) //open the file
    {
      cout << "ERROR: Cannot open the output file: " << projName << endl;
      return;
    }else{
      fwrite(&mcaProjOut,sizeof(mcaProjOut),1,projOut);
      fclose(projOut);
    }

  }

}

void EGamma_PIDsep_mca::SortData(char const *afile, char const *calfile, const unsigned int numP, const unsigned int numA, const int gateEmin, const int gateEmax, double keVPerBin, const int writeProj){

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
  if(AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }
  else
  {
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
      
      if(evtNumProtons==numP){
        if(evtNumAlphas==numA){
          if((evtNumProtons+evtNumAlphas)<=MAX_NUM_PARTICLE){

            for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                if(!add_hit->BGOFired() && add_hit->GetEnergy() > 15){
                  //TIGRESS PID separated addback energy
                  double eGamma1 = add_hit->GetEnergy();
                  if(eGamma1 < gateEmax){
                    if(eGamma1 > gateEmin){
                      for(int tigHitIndAB2=0;tigHitIndAB2<tigress->GetAddbackMultiplicity();tigHitIndAB2++){
                        if(tigHitIndAB != tigHitIndAB2){
                          if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                            add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                            if(!add_hit2->BGOFired() && add_hit2->GetEnergy() > 15){
                              int eGamma2 = (int)(add_hit2->GetEnergy()/keVPerBin);
                              if(eGamma2>=0 && eGamma2<S32K){
                                mcaOut[getTIGRESSRing(add_hit2->GetPosition().Theta()*180./PI)+1][eGamma2]++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  if(writeProj){
                    int projE = (int)(eGamma1/keVPerBin);
                    if(projE>=0 && projE < S32K){
                      mcaProjOut[getTIGRESSRing(add_hit->GetPosition().Theta()*180./PI)+1][projE]++;
                    }
                  }
                  
                }
              }
            }
          }else{
            cout << "Event " << jentry << " has too many charged particles (" << evtNumProtons+evtNumAlphas << ")!" << endl;
          }
        }
      }
        
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  
  analysisfile->Close();
  
}

int main(int argc, char **argv){

  EGamma_PIDsep_mca *mysort = new EGamma_PIDsep_mca();

  char const *afile;
  char const *calfile;
  unsigned int numP, numA;
  int gateEmin, gateEmax;
  double keVPerBin = 1.0;
  int writeProj = 0;
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
  if((argc != 7)&&(argc != 8)&&(argc != 9)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EGamma_PIDsep_mca analysis_tree calibration_file numP numA gateE_min gateE_max keV_per_bin write_proj" << endl;
    cout << "  *analysis_tree* can be a single tree (extension .root) or a" << endl;
    cout << "  plaintext list of trees (one per line, extension .list)." << endl;
    cout << "  *keV_per_bin* defaults to 1." << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else{
    afile = argv[1];
    calfile = argv[2];
    numP = (unsigned int)atoi(argv[3]);
    numA = (unsigned int)atoi(argv[4]);
    gateEmin = atoi(argv[5]);
    gateEmax = atoi(argv[6]);
    if(argc > 7){
      keVPerBin = atof(argv[7]);
      if(argc == 9){
        writeProj = atoi(argv[8]);
      }
    }
  }

  if(gateEmax <= gateEmin){
    cout << "ERROR: Energy gate is zero or negative width!" << endl;
    return 0;
  }
  if((numP+numA)>MAX_NUM_PARTICLE){
    cout << "ERROR: Proton/alpha gate exceeds the maximum possible number of charged particles (" << MAX_NUM_PARTICLE << ")!" << endl;
    return 0;
  }
  if(keVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
    return 0;
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum
  memset(mcaProjOut,0,sizeof(mcaProjOut));

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(strcmp(dot + 1, "root") == 0){

    //single analysis tree
    cout << "Analysis tree: " << afile << endl;
    cout << "Calibration file: " << calfile << endl;
    cout << "Gating on " << numP << " protons, " << numA << " alphas." << endl;
    cout << "Energy gate: [" << gateEmin << " " << gateEmax << "]" << endl;
    cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

    theApp=new TApplication("App", &argc, argv);

    mysort->SortData(afile, calfile, numP, numA, gateEmin, gateEmax, keVPerBin, writeProj);

  }else if(strcmp(dot + 1, "list") == 0){

    //list of analysis trees
    cout << "Analysis tree list: " << afile << endl;
    cout << "Calibration file: " << calfile << endl;
    cout << "Gating on " << numP << " protons, " << numA << " alphas." << endl;
    cout << "Energy gate: [" << gateEmin << " " << gateEmax << "]" << endl;
    cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

    theApp=new TApplication("App", &argc, argv);

    FILE *listfile;
    char str[256];

    if((listfile=fopen(afile,"r"))==NULL){
      cout << "ERROR: Cannot open the input file: " << afile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          mysort->SortData(str, calfile, numP, numA, gateEmin, gateEmax, keVPerBin, writeProj);
        }
      }
    }

  }else{
    cout << "ERROR: Invalid analysis tree file type!" << endl;
    return 0;
  }
  
  mysort->WriteData(numP, numA, gateEmin, gateEmax, writeProj);

  return 0;
}
