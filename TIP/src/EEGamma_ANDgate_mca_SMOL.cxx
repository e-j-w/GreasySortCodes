//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_ANDgate_mca_SMOL_cxx
#include "common.h"
#include "EEGamma_ANDgate_mca_SMOL.h"

using namespace std;

void EEGamma_ANDgate_mca_SMOL::WriteData(const char* outName){

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

void EEGamma_ANDgate_mca_SMOL::SortData(char const *sfile, const double numGates, const double eLow[MAX_GATES], const double eHigh[MAX_GATES], const double keVPerBin){

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

    //check whether AND condition is met
    uint32_t hitpattern = 0;
    uint32_t gateHP = 0;
    uint8_t numGatesMet = 0;
    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){
      if(tigHitIndAB < 32){ //number of bits in hitpattern
        for(int i=0;i<numGates;i++){
          if(i<32){ //number of bits in gate hitpattern
            if(!(gateHP&(1U<<i))){
              //no other hits in this gate 
              if((sortedEvt.tigHit[tigHitIndAB].energy >= eLow[i])&&(sortedEvt.tigHit[tigHitIndAB].energy <= eHigh[i])){
                hitpattern |= (1U<<tigHitIndAB); //flag hit
                gateHP |= (1U << i);
                numGatesMet++;
                break;
              }
            }
          }
        }
      }
    }

    if(numGatesMet == numGates){
      for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){
        if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
          if(!(hitpattern&(1U<<tigHitIndAB))){
            //hit wasn't flagged as being in any gate
            int eGamma = (int)(sortedEvt.tigHit[tigHitIndAB].energy/keVPerBin);
            if(eGamma>=0 && eGamma<S32K){
              double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180./PI;
              mcaOut[getTIGRESSRing(theta)+1][eGamma]++;
              mcaOut[getTIGRESSSegmentRing(theta)+7][eGamma]++;
              mcaOut[0][eGamma]++;
            }
          }
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

  EEGamma_ANDgate_mca_SMOL *mysort = new EEGamma_ANDgate_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow[MAX_GATES], eHigh[MAX_GATES];
  int numGates = 0;
  printf("Starting EEGamma_ANDgate_mca_SMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  if((argc%2 != 0)||(argc < 6)||(argc > (4+(MAX_GATES*2)))){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_ANDgate_mca_SMOL smol_file output_file keV_per_bin EGateLow1 EGateHigh1 EGateLow2 EGateHigh2..." << endl;
    cout << "  Up to " << MAX_GATES << " energy gates may be specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    keVPerBin = atof(argv[3]);
    for(int i=4;i<argc;i+=2){
      if(numGates > MAX_GATES){
        cout << "WARNING: Number of gates exceeds maximum!" << endl;
        break;
      }
      eLow[numGates] = atof(argv[i]);
      eHigh[numGates] = atof(argv[i+1]);
      numGates++;
    }
  }

  if(keVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
    return 0;
  }

  for(int i=0;i<numGates;i++){
    if(eLow[i] > eHigh[i]){
      //swap energy gate bounds
      double tmp = eHigh[i];
      eHigh[i] = eLow[i];
      eLow[i] = tmp;
    }
  }
  

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  for(int i=0;i<numGates;i++){
    cout << "Energy gate " << i << ": [" << eLow[i] << " " << eHigh[i] << "]" << endl;
  }
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

  mysort->SortData(sfile, numGates, eLow, eHigh, keVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
