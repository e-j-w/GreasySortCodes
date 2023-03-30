//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EGamma_mca_SMOL_nopos15_cxx
#include "common.h"
#include "EGamma_mca_SMOL_nopos15.h"

using namespace std;

void EGamma_mca_SMOL_nopos15::WriteData(const char* outName){

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

void EGamma_mca_SMOL_nopos15::SortData(char const *sfile, double keVPerBin){

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
        if((sortedEvt.tigHit[tigHitIndAB].core < 56)||(sortedEvt.tigHit[tigHitIndAB].core >= 60)){
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

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  
  fclose(inp);
  
}

int main(int argc, char **argv){

  EGamma_mca_SMOL_nopos15 *mysort = new EGamma_mca_SMOL_nopos15();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  printf("Starting EGamma_mca_SMOL_nopos15\n");

  if((argc != 3)&&(argc != 4)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EGamma_mca_SMOL_nopos15 smol_file output_file keV_per_bin" << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    if(argc > 3){
      keVPerBin = atof(argv[3]);
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
