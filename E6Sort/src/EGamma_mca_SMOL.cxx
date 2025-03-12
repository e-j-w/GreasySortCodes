//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EGamma_mca_SMOL_cxx
#include "common.cxx"
#include "EGamma_mca_SMOL.h"

using namespace std;

void EGamma_mca_SMOL::WriteData(const char* outName){

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

uint64_t EGamma_mca_SMOL::SortData(char const *sfile, const double keVPerBin){

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

    for(int ABHitInd = 0; ABHitInd < sortedEvt.header.numABHits; ABHitInd++){
      if(sortedEvt.ABHit[ABHitInd].energy > MIN_HPGE_EAB){
        int eGamma = (int)(sortedEvt.ABHit[ABHitInd].energy/keVPerBin);
        if(eGamma>=0 && eGamma<S32K){
          mcaOut[eGamma]++;
        }
      }
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  
  fclose(inp);
  
  return sentries;

}

int main(int argc, char **argv){

  EGamma_mca_SMOL *mysort = new EGamma_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  printf("Starting EGamma_mca_SMOL\n");

  if((argc != 3)&&(argc != 4)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EGamma_mca_SMOL smol_file output_dmca_file keV_per_bin" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
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

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree or list file name." << endl;
    return 0;
  }

  uint64_t numSepEvts = 0U;
  if(strcmp(dot + 1, "smole6") == 0){
    printf("SMOL tree: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, outfile, keVPerBin);
    numSepEvts += mysort->SortData(sfile, keVPerBin);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, outfile, keVPerBin);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, keVPerBin);
        }
      }
    }
  }else{
    cout << "ERROR: improper file extension for SMOL tree or list (should be .smole6 or .list)." << endl;
    return 0;
  }
  
  mysort->WriteData(outfile);
  cout << "Wrote " << numSepEvts << " separated events to: " << outfile << endl;

  return 0;
}
