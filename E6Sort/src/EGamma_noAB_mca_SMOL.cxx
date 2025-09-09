//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EGamma_noAB_mca_SMOL_cxx
#include "common.cxx"
#include "EGamma_noAB_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

void EGamma_noAB_mca_SMOL::WriteData(const char* outName){

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

uint64_t EGamma_noAB_mca_SMOL::SortData(char const *sfile, const double keVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  uint64_t pileupCtrs[16];
  fread(&sentries,sizeof(uint64_t),1,inp);
  uint64_t smolVersion = (uint64_t)(sentries >> 48);
  if(smolVersion > 0){
    fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
    printf("\nNumber of hits of each pileup type:\n");
    uint64_t totalHits = 0;
    for(uint8_t i=0; i<16; i++){
      printf("Pileup type %2u: %Lu\n",i,pileupCtrs[i]);
      totalHits += pileupCtrs[i];
    }
    printf("Total hits:     %Lu\n",totalHits);
    long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
    printf("Fraction of hits with type 1 (no pileup): %Lf\n",frac);
  }
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  sorted_evt sortedEvt;

  //construct 180 degree summing hit map
  memset(hitMap180deg,0,sizeof(hitMap180deg));
  for(uint8_t i=0;i<64;i++){ //first core
    for(uint8_t j=0;j<64;j++){ //coinc core
        if(i!=j){
            if(getGeVector(i,0,1).Angle(getGeVector(j,0,1))*180.0/PI > 175.0){ //same effect for any value down to 165 degrees
                hitMap180deg[i][j] = 1;
                continue; //check the next coinc core
            }
        }
    }
  }

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
        if(eGamma>=0 && eGamma<S32K){
          mcaOut[0][eGamma]++;
        }
      }
      //make summing histograms (for singles)
      for(int noABHitInd3 = noABHitInd+1; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){ //index to avoid double counting in sum histo (B+C 180 sum is the same event as C+B)
        if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd].core & 63U] != 0){
          Double_t tDiffSum = (noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd));
          if((tDiffSum >= SUM_TIMING_GATE_MIN)&&(tDiffSum <= SUM_TIMING_GATE_MAX)){ //timing condition (sum)
            int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
            int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
            int eGammaSum = (int)((sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd3].energy)/keVPerBin);
            //int eGammaSum = (int)(correctSumE(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd3].energy,eLow,eHigh,sumTailLambda,tDiffSum)/keVPerBin);
            if(eGammaSum>=0 && eGammaSum<S32K){
                mcaOut[2][eGammaSum]++; //fill 180 degree sum histogram
            }
            if(eGamma3>=0 && eGamma3<S32K){
                mcaOut[1][eGamma3]++; //fill 180 degree projection histogram
            }
            if(eGamma>=0 && eGamma<S32K){
                mcaOut[1][eGamma]++; //fill 180 degree projection histogram
            }
          }
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

  EGamma_noAB_mca_SMOL *mysort = new EGamma_noAB_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  printf("Starting EGamma_noAB_mca_SMOL\n");

  if((argc != 3)&&(argc != 4)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EGamma_noAB_mca_SMOL smol_file output_dmca_file keV_per_bin" << endl;
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
  if(strcmp(dot + 1, "smol") == 0){
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
    cout << "ERROR: improper file extension for SMOL tree or list (should be .smol or .list)." << endl;
    return 0;
  }
  
  mysort->WriteData(outfile);
  cout << "Wrote " << numSepEvts << " separated events to: " << outfile << endl;

  return 0;
}
