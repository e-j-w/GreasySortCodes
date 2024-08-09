//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_no90deg_mca_SMOL_cxx
#include "common.h"
#include "EEGamma_no90deg_mca_SMOL.h"

using namespace std;

void EEGamma_no90deg_mca_SMOL::WriteData(const char* outName){

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

void EEGamma_no90deg_mca_SMOL::SortData(char const *sfile, const double eLow, const double eHigh, const double keVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){
      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
        if((sortedEvt.tigHit[tigHitIndAB].energy >= eLow)&&(sortedEvt.tigHit[tigHitIndAB].energy <= eHigh)){
          Int_t ring = getTIGRESSRing(getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180./PI);
          if((ring != 2)&&(ring != 3)){
            for(int tigHitIndAB2 = 0; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
              if(tigHitIndAB2 != tigHitIndAB){
                Int_t ring2 = getTIGRESSRing(getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,sortedEvt.tigHit[tigHitIndAB2].seg).Theta()*180./PI);
                if((ring2 != 2)&&(ring2 != 3)){
                  //not at 90 degrees
                  int eGamma = (int)(sortedEvt.tigHit[tigHitIndAB2].energy/keVPerBin);
                  if(eGamma>=0 && eGamma<S32K){
                    double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,sortedEvt.tigHit[tigHitIndAB2].seg).Theta()*180./PI;
                    mcaOut[getTIGRESSRing(theta)+1][eGamma]++;
                    mcaOut[getTIGRESSSegmentRing(theta)+7][eGamma]++;
                    mcaOut[0][eGamma]++;
                  }
                }
              }
            }
          }
          break;
        }
      }
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  
  fclose(inp);
  
}

int main(int argc, char **argv){

  EEGamma_no90deg_mca_SMOL *mysort = new EEGamma_no90deg_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  printf("Starting EEGamma_no90deg_mca_SMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  if((argc != 5)&&(argc != 6)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_no90deg_mca_SMOL smol_file EGateLow EGateHigh output_file keV_per_bin" << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    eLow = atof(argv[2]);
    eHigh = atof(argv[3]);
    outfile = argv[4];
    if(argc > 5){
      keVPerBin = atof(argv[5]);
    }
  }

  if(keVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
    return 0;
  }

  if(eLow > eHigh){
    //swap energy gate bounds
    double tmp = eHigh;
    eHigh = eLow;
    eLow = tmp;
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  cout << "Energy gate: [" << eLow << " " << eHigh << "]" << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

  mysort->SortData(sfile, eLow, eHigh, keVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
