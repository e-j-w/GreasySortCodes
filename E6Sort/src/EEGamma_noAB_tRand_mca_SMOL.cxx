//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_noAB_tRand_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_noAB_tRand_mca_SMOL.h"

using namespace std;

void EEGamma_noAB_tRand_mca_SMOL::WriteData(const char* outName){

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

void EEGamma_noAB_tRand_mca_SMOL::SortData(char const *sfile, const double eLow, const double eHigh, const double keVPerBin){

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

    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        if((sortedEvt.noABHit[noABHitInd].energy >= eLow)&&(sortedEvt.noABHit[noABHitInd].energy <= eHigh)){
          for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            //perform time separation
            Double_t tDiff = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2);
            if((tDiff >= hpgehpgeTRandGate[0])&&(tDiff <= hpgehpgeTRandGate[1])){
              int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
              if(eGamma>=0 && eGamma<S32K){
                mcaOut[eGamma]++;
              }
            }
          }
          //shouldn't break, what if there are 2 gammas in the gate?
          //break;
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

  EEGamma_noAB_tRand_mca_SMOL *mysort = new EEGamma_noAB_tRand_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  printf("Starting EEGamma_noAB_tRand_mca_SMOL\n");

  if((argc != 5)&&(argc != 6)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_noAB_tRand_mca_SMOL smol_file EGateLow EGateHigh output_file keV_per_bin" << endl;
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
