//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_noAB_mca_SMOL_splitruns_cxx
#include "common.cxx"
#include "EEGamma_noAB_mca_SMOL_splitruns.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint64_t filehits;
uint32_t numFilesWritten;
uint64_t hitsFilled[9]; //bit-pattern desrcibing the hits filled in each spectrum

void WriteData(const char* outName){

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

uint64_t EEGamma_noAB_mca_SMOL_splitruns::SortData(const char *sfile, const char *outfile, const double eLow, const double eHigh, const double keVPerBin, const uint8_t discardPileup, const uint64_t entriesPerFile){

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
  srand(92370104);

  //allow decimation of sorted events (for debugging/tuning)
  Long64_t increment = 1;
  if(increment == 1){
    printf("\nSorting events...\n");
  }else if(increment > 0){
    printf("\nSorting every %i events...\n",increment);
  }else{
    increment = 1;
  }

  for(Long64_t jentry = 0; jentry < sentries; jentry+=increment){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    if(sortedEvt.header.numNoABHits <= 64){

      //hits can fit in uint64_t bit-pattern
      memset(hitsFilled,0,sizeof(hitsFilled));

      for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

        if(discardPileup && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
          continue; //skip pileup hit
        }

        //fill singles spectra
        //int singlesE = (int)( (timeRandomOffsetFactor*sortedEvt.noABHit[noABHitInd].energy*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
        int singlesE = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
        if((singlesE>=0) && (singlesE<S32K)){
          mcaOut[6][singlesE]++;
        }
        //make summing histograms (for singles)
        for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
          if(noABHitInd3 != noABHitInd){
            if(discardPileup && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
              continue; //skip pileup hit
            }
            if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd].core & 63U] != 0){
              Double_t tDiffSum = noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd);
              if((tDiffSum >= SUM_TIMING_GATE_MIN)&&(tDiffSum <= SUM_TIMING_GATE_MAX)){ //timing condition (sum)
                int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
                int eGammaSum = (int)((sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd3].energy)/keVPerBin);
                //int eGammaSum = (int)(correctSumE(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd3].energy,tDiffSum)/keVPerBin);
                if(eGammaSum>=0 && eGammaSum<S32K){
                  if(!(hitsFilled[8] & ((uint64_t)(1) << noABHitInd))){
                    if(!(hitsFilled[8] & ((uint64_t)(1) << noABHitInd3))){
                      hitsFilled[8] |= ((uint64_t)(1) << noABHitInd);
                      hitsFilled[8] |= ((uint64_t)(1) << noABHitInd3);
                      mcaOut[8][eGammaSum]++; //fill 180 degree sum histogram
                    }
                  }
                }
                if(eGamma3>=0 && eGamma3<S32K){
                  if(!(hitsFilled[7] & ((uint64_t)(1) << noABHitInd3))){
                    hitsFilled[7] |= ((uint64_t)(1) << noABHitInd3);
                    mcaOut[7][eGamma3]++; //fill 180 degree projection histogram
                  }
                }
                if(eGamma>=0 && eGamma<S32K){
                  if(!(hitsFilled[7] & ((uint64_t)(1) << noABHitInd))){
                    hitsFilled[7] |= ((uint64_t)(1) << noABHitInd);
                    mcaOut[7][eGamma]++; //fill 180 degree projection histogram
                  }
                }
              }
            }
          }
        }
        
        if((sortedEvt.noABHit[noABHitInd].energy >= eLow)&&(sortedEvt.noABHit[noABHitInd].energy <= eHigh)){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            if(discardPileup && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
              continue; //skip pileup hit
            }
            //if(noABHitInd != noABHitInd2){
            if(((sortedEvt.noABHit[noABHitInd].core & 63U)/4)!=((sortedEvt.noABHit[noABHitInd2].core & 63U)/4)){ //try to reduce crosstalk... doesn't seem to do anything regarding sum peak shapes, but seems to align time-random with singles data
              Double_t tDiff = (noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
              uint8_t numCFDFail = 0;
              if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(((numCFDFail % 2 == 0)&&(tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX))||((numCFDFail % 2 == 0)&&(tDiff >= COINC_TIMING_GATE_CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_CFDFAIL_MAX))){
                if(!(hitsFilled[0] & ((uint64_t)(1) << noABHitInd2))){
                  int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                  if(eGamma>=0 && eGamma<S32K){
                    hitsFilled[0] |= ((uint64_t)(1) << noABHitInd2);
                    mcaOut[0][eGamma]++; //fill true coincidence histogram
                  }
                }
              }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                //time random
                if(!(hitsFilled[3] & ((uint64_t)(1) << noABHitInd2))){
                  int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                  //double eGamma = ( (timeRandomOffsetFactor*sortedEvt.noABHit[noABHitInd2].energy*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                  if(((int)eGamma)>=0 && ((int)eGamma)<S32K){
                    hitsFilled[3] |= ((uint64_t)(1) << noABHitInd2);
                    mcaOut[3][((int)eGamma)]++; //fill time-random 'coincidence' histogram
                  }
                }
              }
              //make summing histograms
              for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                if(discardPileup && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
                  continue; //skip pileup hit
                }
                //if(noABHitInd3 != noABHitInd){
                if(((sortedEvt.noABHit[noABHitInd].core & 63U)/4)!=((sortedEvt.noABHit[noABHitInd3].core & 63U)/4)){
                  if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                    Double_t tDiffSum = (noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd2));
                    if((tDiffSum >= SUM_TIMING_GATE_MIN)&&(tDiffSum <= SUM_TIMING_GATE_MAX)){ //timing condition (sum)
                      Double_t tDiffFirstSumHit = (noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
                      uint8_t numCFDFail = 0;
                      if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                        numCFDFail++;
                      }
                      if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                        numCFDFail++;
                      }
                      if(((numCFDFail % 2 == 0)&&(tDiffFirstSumHit >= COINC_TIMING_GATE_MIN)&&(tDiffFirstSumHit <= COINC_TIMING_GATE_MAX))||((numCFDFail % 2 == 1)&&(tDiffFirstSumHit >= COINC_TIMING_GATE_CFDFAIL_MIN)&&(tDiffFirstSumHit <= COINC_TIMING_GATE_CFDFAIL_MAX))){ //timing condition (original energy gate)
                        //printf("tDiffSum: %f\n",tDiffSum);
                        int eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                        int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
                        int eGammaSum = (int)((sortedEvt.noABHit[noABHitInd2].energy + sortedEvt.noABHit[noABHitInd3].energy)/keVPerBin);
                        //int eGammaSum = (int)(correctSumE(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd3].energy,tDiffSum)/keVPerBin);
                        if(eGammaSum>=0 && eGammaSum<S32K){
                          if(!(hitsFilled[2] & ((uint64_t)(1) << noABHitInd2))){
                            if(!(hitsFilled[2] & ((uint64_t)(1) << noABHitInd3))){
                              hitsFilled[2] |= ((uint64_t)(1) << noABHitInd2);
                              hitsFilled[2] |= ((uint64_t)(1) << noABHitInd3);
                              mcaOut[2][eGammaSum]++; //fill 180 degree sum histogram
                            }
                          }
                        }
                        if(eGamma3>=0 && eGamma3<S32K){
                          if(!(hitsFilled[1] & ((uint64_t)(1) << noABHitInd3))){
                            hitsFilled[1] |= ((uint64_t)(1) << noABHitInd3);
                            mcaOut[1][eGamma3]++; //fill 180 degree projection histogram
                          }
                        }
                        if(eGamma>=0 && eGamma<S32K){
                          if(!(hitsFilled[1] & ((uint64_t)(1) << noABHitInd2))){
                            hitsFilled[1] |= ((uint64_t)(1) << noABHitInd2);
                            mcaOut[1][eGamma]++; //fill 180 degree projection histogram
                          }
                        }
                      }else if((tDiffFirstSumHit >= TRANDOM_GATE_MIN)&&(tDiffFirstSumHit <= TRANDOM_GATE_MAX)){
                        double eGamma = ( (sortedEvt.noABHit[noABHitInd2].energy)/keVPerBin );
                        double eGamma3 = ( (sortedEvt.noABHit[noABHitInd3].energy)/keVPerBin );
                        int eGammaSum = (int)((sortedEvt.noABHit[noABHitInd2].energy + sortedEvt.noABHit[noABHitInd3].energy)/keVPerBin);
                        //double eGamma = ( (timeRandomOffsetFactor*sortedEvt.noABHit[noABHitInd2].energy*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                        //double eGamma3 = ( (timeRandomOffsetFactor*sortedEvt.noABHit[noABHitInd3].energy*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                        //double eGammaSum = (correctSumE(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd3].energy,tDiffSum)*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin;
                        if(((int)eGammaSum)>=0 && ((int)eGammaSum)<S32K){
                          if(!(hitsFilled[5] & ((uint64_t)(1) << noABHitInd2))){
                            if(!(hitsFilled[5] & ((uint64_t)(1) << noABHitInd3))){
                              hitsFilled[5] |= ((uint64_t)(1) << noABHitInd2);
                              hitsFilled[5] |= ((uint64_t)(1) << noABHitInd3);
                              mcaOut[5][((int)eGammaSum)]++; //fill 180 degree sum histogram (time random)
                            }
                          }
                        }
                        if(((int)eGamma3)>=0 && ((int)eGamma3)<S32K){
                          if(!(hitsFilled[4] & ((uint64_t)(1) << noABHitInd3))){
                            hitsFilled[4] |= ((uint64_t)(1) << noABHitInd3);
                            mcaOut[4][((int)eGamma3)]++; //fill 180 degree projection histogram (time random)
                          }
                        }
                        if(((int)eGamma)>=0 && ((int)eGamma)<S32K){
                          if(!(hitsFilled[4] & ((uint64_t)(1) << noABHitInd2))){
                            hitsFilled[4] |= ((uint64_t)(1) << noABHitInd2);
                            mcaOut[4][((int)eGamma)]++; //fill 180 degree projection histogram (time random)
                          }
                        }
                      }
                    }
                  }
                }
              }
              //}
            }
          }
          //shouldn't break, what if there are 2 gammas in the gate?
          //break;
          
        }
      }

      filehits+=sortedEvt.header.numNoABHits;

    }

    if (jentry % 90713 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;

    if(filehits >= entriesPerFile){
        char fileName[256];
        snprintf(fileName,256,"%s_%u.dmca",outfile,numFilesWritten);
        cout << endl << "Writing " << filehits << " hits to file: " << fileName << endl << endl;
        WriteData(fileName); //write data to disk
        memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum
        filehits = 0; //zero out the number of entries sorted for the next file
        numFilesWritten++;
    }
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  
  fclose(inp);

  return sentries;
  
}

int main(int argc, char **argv){

  EEGamma_noAB_mca_SMOL_splitruns *mysort = new EEGamma_noAB_mca_SMOL_splitruns();

  const char *sfile;
  const char *outfile;
  uint8_t discardPileup = 0;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  uint64_t entriesPerFile = 0;
  printf("Starting EEGamma_noAB_mca_SMOL_splitruns\n");

  if((argc != 6)&&(argc != 7)&&(argc != 8)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_noAB_mca_SMOL_splitruns smol_file_list EGateLow EGateHigh output_dmca_file_prefix events_per_split keV_per_bin discard_pileup" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    cout << "  *discard_pileup* can be either 0 (false, default if not specified) or 1 (true)." << endl;
    return 0;
  }else{
    sfile = argv[1];
    eLow = atof(argv[2]);
    eHigh = atof(argv[3]);
    outfile = argv[4];
    entriesPerFile = (uint64_t)atoll(argv[5]);
    if(argc > 6){
      keVPerBin = atof(argv[6]);
      if(argc > 7){
        if(atoi(argv[7]) != 0){
          discardPileup = 1;
        }
      }

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
  filehits = 0; //zero out the number of entries sorted for the next file
  numFilesWritten = 0;

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree list file name." << endl;
    return 0;
  }

  uint64_t numSepEvts = 0U;
  if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nEnergy gate: [%0.2f %0.2f]\nOutput file: %s\nEntries per file: %lu\n%0.2f keV per bin\n", sfile, eLow, eHigh, outfile, entriesPerFile, keVPerBin);
    if(discardPileup == 1){
      printf("Discarding pileup hits.\n");
    }

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
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, outfile, eLow, eHigh, keVPerBin, discardPileup, entriesPerFile);
        }
      }
    }
  }else{
    cout << "ERROR: improper file extension for SMOL tree list (should be .list)." << endl;
    return 0;
  }
  
  char fileName[256];
  snprintf(fileName,256,"%s_%u.dmca",outfile,numFilesWritten);
  cout << "Writing " << filehits << " events to file: " << fileName << endl;
  WriteData(fileName); //write remaining data to disk

  cout << "Sorted a total of " << numSepEvts << " events." << endl;

  return 0;
}
