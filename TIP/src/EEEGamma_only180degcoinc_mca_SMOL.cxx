//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEEGamma_only180degcoinc_mca_SMOL_cxx
#include "common.h"
#include "EEEGamma_only180degcoinc_mca_SMOL.h"

using namespace std;

void EEEGamma_only180degcoinc_mca_SMOL::WriteData(const char* outName){

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

void EEEGamma_only180degcoinc_mca_SMOL::SortData(char const *sfile, const double eLow, const double eHigh, const double keVPerBin){

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

    if(sortedEvt.header.numTigHits >= 3){
        for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){

        int eGamma = (int)(sortedEvt.tigHit[tigHitIndAB].energy/keVPerBin);
        if((eGamma >= eLow)&&(eGamma <= eHigh)){
            //TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,0);
            for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
                for(int tigHitIndAB3 = tigHitIndAB2+1; tigHitIndAB3 < sortedEvt.header.numTigHits; tigHitIndAB3++){
                    TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,0);
                    TVector3 vec3 = getTigVector(sortedEvt.tigHit[tigHitIndAB3].core,0);
                    Double_t angle = vec2.Angle(vec3)*180./PI; //angle between hits
                    if((angle >= 170)&&(angle <= 190)){
                        //printf("Angle: %f\n",angle);
                        //coincidence at 180 degrees
                        int eGamma2 = (int)(sortedEvt.tigHit[tigHitIndAB2].energy/keVPerBin);
                        if(eGamma2>=0 && eGamma2<S32K){
                            double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,sortedEvt.tigHit[tigHitIndAB2].seg).Theta()*180./PI;
                            mcaOut[getTIGRESSRing(theta)+1][eGamma2]++;
                            mcaOut[getTIGRESSSegmentRing(theta)+7][eGamma2]++;
                            mcaOut[0][eGamma2]++;
                        }
                        int eGamma3 = (int)(sortedEvt.tigHit[tigHitIndAB3].energy/keVPerBin);
                        if(eGamma3>=0 && eGamma3<S32K){
                            double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB3].core,sortedEvt.tigHit[tigHitIndAB3].seg).Theta()*180./PI;
                            mcaOut[getTIGRESSRing(theta)+1][eGamma3]++;
                            mcaOut[getTIGRESSSegmentRing(theta)+7][eGamma3]++;
                            mcaOut[0][eGamma3]++;
                        }
                        break;
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

  EEEGamma_only180degcoinc_mca_SMOL *mysort = new EEEGamma_only180degcoinc_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  printf("Starting EEEGamma_only180degcoinc_mca_SMOL\n");
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
    cout << "User provides a gamma gate, spectra will contain only events in which there are 2 coincident gammas, and those two are emitted at 180 degrees (back-to-back). This is used for summing-out corrections." << endl;
    cout << "Arguments: EEEGamma_only180degcoinc_mca_SMOL smol_file EGateLow EGateHigh output_file keV_per_bin" << endl;
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
