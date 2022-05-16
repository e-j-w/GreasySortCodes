//Plots TIGRESS Doppler corrected gamma-gamma matrices for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEDopp_PIDsep_mca_cxx
#include "common.h"
#include "EDopp_PIDsep_mca.h"

using namespace std;

void EDopp_PIDsep_mca::SortData(char const *afile, char const *calfile, char const *outfile, const unsigned int numP, const unsigned int numA, const int gateEmin, const int gateEmax){

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

  double_t tipPID = -1000.0;
  unsigned int evtNumProtons, evtNumAlphas;
  bool suppAdd = false;

  //Defining Pointers
  TTipHit *tip_hit;
  TTigressHit *add_hit, *add_hit2;
  const std::vector<Short_t> *wf; //for CsI waveform

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  printf("\nSorting analysis events...\n");
  for (int jentry = 0; jentry < analentries; jentry++){

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

      uint32_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup

      //count the number of protons or alphas
      for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
        if(passedtimeGate&(1U<<tipHitInd)){
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
              if(passedtimeGate&(1U<<(tigHitIndAB+16))){
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                suppAdd = add_hit->BGOFired();
                //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                if(!suppAdd && add_hit->GetEnergy() > 15){
                  //TIGRESS PID separated addback energy
                  double eDopp1 = getEDoppFusEvap(add_hit,tip,passedtimeGate,gates);
                  if(eDopp1 < gateEmax){
                    if(eDopp1 > gateEmin){
                      for(int tigHitIndAB2=tigHitIndAB+1;tigHitIndAB2<tigress->GetAddbackMultiplicity();tigHitIndAB2++){
                        add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                        suppAdd = add_hit->BGOFired();
                        if(!suppAdd && add_hit2->GetEnergy() > 15){
                          int eDopp2 = (int)(getEDoppFusEvap(add_hit2,tip,passedtimeGate,gates)*2);
                          mcaOut[getTIGRESSRing(add_hit->GetPosition().Theta()*180./PI)][eDopp2]++;
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
        }
      }
        
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  FILE *out;
  if((out = fopen(outfile, "w")) == NULL) //open the file
  {
    cout << "ERROR: Cannot open the output file: " << outfile << endl;
    return;
  }else{
    fwrite(&mcaOut,sizeof(mcaOut),1,out);
    fclose(out);
  }
  
  
}

int main(int argc, char **argv){

  EDopp_PIDsep_mca *mysort = new EDopp_PIDsep_mca();

  char const *afile;
  char const *outfile;
  char const *calfile;
  unsigned int numP, numA;
  int gateEmin, gateEmax;
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
  if(argc != 8){
    cout << "Generates TIGRESS mca spectra for Doppler corrected, PID and time separated data." << endl;
    cout << "Arguments: EEDopp_PIDsep_mca analysis_tree calibration_file numP numA gateE_min gateE_max output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else if(argc == 4){
    afile = argv[1];
    calfile = argv[2];
    numP = (unsigned int)atoi(argv[3]);
    numA = (unsigned int)atoi(argv[4]);
    gateEmin = atoi(argv[5]);
    gateEmax = atoi(argv[6]);
    outfile = argv[7];
  }

  if(gateEmax <= gateEmin){
    cout << "ERROR: Energy gate is zero or negative width!" << endl;
    return 0;
  }
  if((numP+numA)>MAX_NUM_PARTICLE){
    cout << "ERROR: Proton/alpha gate exceeds the maximum possible number of charged particles (" << MAX_NUM_PARTICLE << ")!" << endl;
    return 0;
  }

  cout << "Analysis file: " << afile << endl;
  cout << "Calibration file: " << calfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Gating on " << numP << " protons, " << numA << " alphas." << endl;
  cout << "Energy gate: [" << gateEmin << " " << gateEmax << "]" << endl;

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile, numP, numA, gateEmin, gateEmax);

  return 0;
}
