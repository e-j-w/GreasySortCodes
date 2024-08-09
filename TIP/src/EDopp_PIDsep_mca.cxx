//Generates TIGRESS Doppler corrected gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EDopp_PIDsep_mca_cxx
#include "common.h"
#include "EDopp_PIDsep_mca.h"

using namespace std;

void EDopp_PIDsep_mca::WriteData(const unsigned int numP, const unsigned int numA){

  char const *outName = Form("EGamma_%up%ua.fmca",numP,numA);
  cout << "Writing gated histogram to: " << outName << endl;

  FILE *out;
  if((out = fopen(outName, "w")) == NULL){ //open the file
    cout << "ERROR: Cannot open the output file: " << outName << endl;
    return;
  }else{
    fwrite(&mcaOut,sizeof(mcaOut),1,out);
    fclose(out);
  }

}

void EDopp_PIDsep_mca::SortData(char const *afile, char const *calfile, const unsigned int numP, const unsigned int numA, double keVPerBin){

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
  TTigressHit *add_hit;

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

    if(tigress){

      if(!tip && ((numP+numA)>0)){
        cout << "WARNING: event " << jentry << " has no TIP data!" << endl;
        continue;
      }
      if(tip && tip->GetMultiplicity()>MAXNUMTIPHIT){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }
      if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }

      uint64_t passedtimeGate = passesTimeGate(tigress,tip,1,numP+numA); //also rejects pileup

      //count the number of protons or alphas
      if((numP+numA)>0){
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
      }
      
      if(evtNumProtons==numP){
        if(evtNumAlphas==numA){
          if((evtNumProtons+evtNumAlphas)<=MAX_NUM_PARTICLE){

            for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                if(!add_hit->BGOFired() && add_hit->GetEnergy() > MIN_TIG_EAB){
                  int eDopp = (int)getEDoppFusEvap(add_hit,tip,passedtimeGate,gates);
                  if(eDopp>=0 && eDopp<S32K){
                    double_t thetaDeg = add_hit->GetPosition().Theta()*180./PI;
                    mcaOut[getTIGRESSRing(thetaDeg)+1][eDopp]++;
                    mcaOut[getTIGRESSSegmentRing(thetaDeg)+7][eDopp]++;
                    mcaOut[0][eDopp]++;
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

  EDopp_PIDsep_mca *mysort = new EDopp_PIDsep_mca();

  char const *afile;
  char const *calfile;
  unsigned int numP, numA;
  double keVPerBin = 1.0;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  if((argc != 5)&&(argc != 6)){
    cout << "Generates TIGRESS Doppler corrected mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EDopp_PIDsep_mca analysis_tree calibration_file numP numA keV_per_bin" << endl;
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
    if(argc > 5){
      keVPerBin = atof(argv[5]);
    }
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

  //Setup TIP PID gates
  cout << "Creating PID gates." << endl;
  gates = new PIDGates;

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(strcmp(dot + 1, "root") == 0){

    //single analysis tree
    cout << "Analysis tree: " << afile << endl;
    cout << "Calibration file: " << calfile << endl;
    cout << "Gating on " << numP << " protons, " << numA << " alphas." << endl;
    cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

    theApp=new TApplication("App", &argc, argv);

    mysort->SortData(afile, calfile, numP, numA, keVPerBin);

  }else if(strcmp(dot + 1, "list") == 0){

    //list of analysis trees
    cout << "Analysis tree list: " << afile << endl;
    cout << "Calibration file: " << calfile << endl;
    cout << "Gating on " << numP << " protons, " << numA << " alphas." << endl;
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
          mysort->SortData(str, calfile, numP, numA, keVPerBin);
        }
      }
    }

  }else{
    cout << "ERROR: Invalid analysis tree file type!" << endl;
    return 0;
  }
  
  mysort->WriteData(numP, numA);

  return 0;
}
