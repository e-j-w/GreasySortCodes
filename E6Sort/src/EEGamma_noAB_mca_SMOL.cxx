//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_noAB_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_noAB_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

//make sure srand() is called somewhere
double rand_expo(double lambda){
  double u;
  u = rand() / (RAND_MAX + 1.0);
  return -log(1- u) / lambda;
}

//function attempting to correct the sum peak energy for 180 degree coincidences
//in true summing, the energy will be shifted down slightly, since the first
//pulse will decay partially before the second hit occurs
float correctSumE(float en, double tDiffSum, double gain, double offset){
  double corrVal = rand_expo(0.50)*tDiffSum/480.0; //divisor more significant for 1740 line, expo more significant for 2950 line
  //printf("corrVal: %f\n",corrVal);
  return (float)((en - corrVal)*gain + offset);
  //lambda = 1.0: not enough skew
  //lambda = 0.5: closer but still not enough skew
  //lambda = 0.25: good for coincident summing, too fast for time random
  //0.35 - slightly too slow for time-random
}

void EEGamma_noAB_mca_SMOL::WriteData(const char* outName){

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



uint64_t EEGamma_noAB_mca_SMOL::SortData(const char *sfile, const double eLow, const double eHigh, const double keVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  srand(92370104);

  //construct 180 degree summing hit map
  memset(hitMap180deg,0,sizeof(hitMap180deg));
  for(uint8_t i=0;i<64;i++){ //first core
      for(uint8_t j=0;j<64;j++){ //coinc core
          if(i!=j){
              if(getGeVector(i,0,1).Angle(getGeVector(j,0,1))*180.0/PI > 175.0){
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
      if((sortedEvt.noABHit[noABHitInd].energy >= eLow)&&(sortedEvt.noABHit[noABHitInd].energy <= eHigh)){
        for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
          if(noABHitInd != noABHitInd2){
            Double_t tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
            if(tDiff <= COINC_TIMING_GATE){
              int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
              if(eGamma>=0 && eGamma<S32K){
                mcaOut[0][eGamma]++; //fill true coincidence histogram
              }
            }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
              //time random
              int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
              if(eGamma>=0 && eGamma<S32K){
                mcaOut[3][eGamma]++; //fill time-random 'coincidence' histogram
              }
            }
            //make summing histograms
            for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
              if((noABHitInd3 != noABHitInd)&&((noABHitInd3 != noABHitInd2))){
                if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                  Double_t tDiffSum = fabs(noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd2));
                  if(tDiffSum <= SUM_TIMING_GATE){ //timing condition (sum)
                    //for comparing the summed hits to the original gated hit, need to consider the timing of the first
                    //sum hit, since in the actual sum peak, the trigger is based off the initial rise of the pulse, 
                    //which will occur at the time of the first hit
                    Double_t tDiffFirstSumHit = tDiff;
                    if(noABHitTime(&sortedEvt,noABHitInd3) < noABHitTime(&sortedEvt,noABHitInd2)){
                      tDiffFirstSumHit = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd3));
                    }
                    if(tDiffFirstSumHit <= COINC_TIMING_GATE){ //timing condition (original energy gate)
                      int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                      int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
                      //int eGammaSum = eGamma + eGamma3;
                      int eGammaSum = (int)(correctSumE(sortedEvt.noABHit[noABHitInd2].energy + sortedEvt.noABHit[noABHitInd3].energy,tDiffSum,1.0,0.0)/keVPerBin);
                      if(eGammaSum>=0 && eGammaSum<S32K){
                          mcaOut[2][eGammaSum]++; //fill 180 degree sum histogram
                      }
                      if(eGamma3>=0 && eGamma3<S32K){
                          mcaOut[1][eGamma3]++; //fill 180 degree projection histogram
                      }
                      if(eGamma>=0 && eGamma<S32K){
                          mcaOut[1][eGamma]++; //fill 180 degree projection histogram
                      }
                    }else if((tDiffFirstSumHit >= TRANDOM_GATE_MIN)&&(tDiffFirstSumHit <= TRANDOM_GATE_MAX)){
                      int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                      int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
                      //int eGammaSum = eGamma + eGamma3;
                      int eGammaSum = (int)(correctSumE(sortedEvt.noABHit[noABHitInd2].energy + sortedEvt.noABHit[noABHitInd3].energy,tDiffSum,1.0,0.0)/keVPerBin);
                      if(eGammaSum>=0 && eGammaSum<S32K){
                          mcaOut[5][eGammaSum]++; //fill 180 degree sum histogram (time random)
                      }
                      if(eGamma3>=0 && eGamma3<S32K){
                          mcaOut[4][eGamma3]++; //fill 180 degree projection histogram (time random)
                      }
                      if(eGamma>=0 && eGamma<S32K){
                          mcaOut[4][eGamma]++; //fill 180 degree projection histogram (time random)
                      }
                    }
                  }
                }
              }
            }
          }
        }
        //shouldn't break, what if there are 2 gammas in the gate?
        //break;
        
      }
    }

    if (jentry % 9713 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  
  fclose(inp);

  return sentries;
  
}

int main(int argc, char **argv){

  EEGamma_noAB_mca_SMOL *mysort = new EEGamma_noAB_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  printf("Starting EEGamma_noAB_mca_SMOL\n");

  if((argc != 5)&&(argc != 6)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_noAB_mca_SMOL smol_file EGateLow EGateHigh output_dmca_file keV_per_bin" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
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

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree or list file name." << endl;
    return 0;
  }

  uint64_t numSepEvts = 0U;
  if(strcmp(dot + 1, "smole6") == 0){
    printf("SMOL tree: %s\nEnergy gate: [%0.2f %0.2f]\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, outfile, keVPerBin);
    numSepEvts += mysort->SortData(sfile, eLow, eHigh, keVPerBin);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nEnergy gate: [%0.2f %0.2f]\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, outfile, keVPerBin);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, eLow, eHigh, keVPerBin);
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
