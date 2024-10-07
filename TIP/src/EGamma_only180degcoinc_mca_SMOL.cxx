//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h

#define EGamma_only180degcoinc_mca_SMOL_cxx
#include "common.h"
#include "EGamma_only180degcoinc_mca_SMOL.h"

using namespace std;

void EGamma_only180degcoinc_mca_SMOL::WriteData(const char* outName){

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

void EGamma_only180degcoinc_mca_SMOL::SortData(char const *sfile, const double keVPerBin){

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

      int eGamma = (int)(sortedEvt.tigHit[tigHitIndAB].energy/keVPerBin);
      if(eGamma>=0 && eGamma<S32K){
        TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,0);
        for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
          TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,0);
          Double_t angle = vec1.Angle(vec2)*180./PI; //angle between hits
          if((angle >= 170)&&(angle <= 190)){
            //printf("Angle: %f\n",angle);
            //coincidence at 180 degrees
            int eGamma2 = (int)(sortedEvt.tigHit[tigHitIndAB2].energy/keVPerBin);
            if(eGamma2>=0 && eGamma2<S32K){
              double theta1 = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180./PI;
              double theta2 = getTigVector(sortedEvt.tigHit[tigHitIndAB2].core,sortedEvt.tigHit[tigHitIndAB2].seg).Theta()*180./PI;
              mcaOut[getTIGRESSRing(theta1)+1][eGamma]++;
              mcaOut[getTIGRESSSegmentRing(theta1)+7][eGamma]++;
              mcaOut[0][eGamma]++;
              mcaOut[getTIGRESSRing(theta2)+1][eGamma]++;
              mcaOut[getTIGRESSSegmentRing(theta2)+7][eGamma]++;
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

  EGamma_only180degcoinc_mca_SMOL *mysort = new EGamma_only180degcoinc_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  printf("Starting EGamma_only180degcoinc_mca_SMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  if((argc != 3)&&(argc != 4)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Only gammas in 180-degree coincidence with another gamma will be included." << endl;
    cout << "Arguments: EGamma_only180degcoinc_mca_SMOL smol_file output_file keV_per_bin" << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    if(argc > 3){
      keVPerBin = atof(argv[4]);
    }
  }

  if(keVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
    return 0;
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

  mysort->SortData(sfile, keVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
