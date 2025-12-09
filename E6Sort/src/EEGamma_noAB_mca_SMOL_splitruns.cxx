//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_noAB_mca_SMOL_splitruns_cxx
#include "common.cxx"
#include "EEGamma_noAB_mca_SMOL_splitruns.h"

using namespace std;

enum sp_enum{
SP_GATED, SP_SUMOUT, SP_SUMIN,
SP_TR_GATED, SP_TR_SUMOUT, SP_TR_SUMIN,
SP_SINGLES, SP_SINGLES_SUMOUT, SP_SINGLES_SUMIN,
SP_ENUM_LENGTH
};

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs
uint64_t filehits;
uint32_t numFilesWritten;



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

uint64_t EEGamma_noAB_mca_SMOL_splitruns::SortData(const char *sfile, const char *outfile, const double eLow, const double eHigh, const double keVPerBin, const uint8_t discardPileup, const uint64_t entriesPerFile, const double offset, const double gain, const double quad){

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

      for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
          continue; //skip pileup hit
        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
          continue; //skip non-pileup hit
        }

        const double hit1E = offset + sortedEvt.noABHit[noABHitInd].energy*gain + sortedEvt.noABHit[noABHitInd].energy*sortedEvt.noABHit[noABHitInd].energy*quad;

        //fill singles spectra
        //int singlesE = (int)( (timeRandomOffsetFactor*hit1E*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
        int singlesE = (int)(hit1E/keVPerBin);
        if((singlesE>=0) && (singlesE<S32K)){
          mcaOut[SP_SINGLES][singlesE]++;
          //filehits++;
        }
        //make summing histograms (for singles)
        for(int noABHitInd3 = noABHitInd+1; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
          if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
            continue; //skip pileup hit
          }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7)))){
            continue; //skip non-pileup hit
          }
          if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd].core & 63U] != 0){
            Double_t tDiffSum = fabs(sortedEvt.noABHit[noABHitInd3].tsDiff - sortedEvt.noABHit[noABHitInd].tsDiff)*10.0;
            //printf("tDiffSum: %f\n",tDiffSum);
            if((tDiffSum >= SUM_TIMING_GATE_MIN)&&(tDiffSum <= SUM_TIMING_GATE_MAX)){ //timing condition (sum)
              //printf("Passed timing condition.\n");

              const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;

              int eGamma3 = (int)(hit3E/keVPerBin);
              int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd].energy*sortedEvt.noABHit[noABHitInd].energy)+(sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
              //int eGammaSum = (int)(correctSumE(hit1E,hit3E,tDiffSum)/keVPerBin);
              if(eGammaSum>=0 && eGammaSum<S32K){
                mcaOut[SP_SINGLES_SUMIN][eGammaSum]++; //fill 180 degree sum histogram
              }
              if(eGamma3>=0 && eGamma3<S32K){
                mcaOut[SP_SINGLES_SUMOUT][eGamma3]++; //fill 180 degree projection histogram
              }
              if(singlesE>=0 && singlesE<S32K){
                mcaOut[SP_SINGLES_SUMOUT][singlesE]++; //fill 180 degree projection histogram
              }
            }
          }
        }

        //custom hit counting condition (1477+685 coinc)
        if((hit1E >= 1466.0)&&(hit1E <= 1486.0)){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
              continue; //skip pileup hit
            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
              continue; //skip non-pileup hit
            }
            const double hit2E = offset + sortedEvt.noABHit[noABHitInd2].energy*gain + sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy*quad;
            if((hit2E >= 674.0)&&(hit2E <= 696.0)){
              //if(((sortedEvt.noABHit[noABHitInd].core & 63U)/4)!=((sortedEvt.noABHit[noABHitInd2].core & 63U)/4)){ //try to reduce crosstalk... doesn't seem to do anything regarding sum peak shapes, but seems to align time-random with singles data
              Double_t tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
              uint8_t numCFDFail = 0;
              if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(((numCFDFail == 0)&&(tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX))||((numCFDFail == 1)&&(tDiff >= COINC_TIMING_GATE_1CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_1CFDFAIL_MAX))||((numCFDFail == 2)&&(tDiff >= COINC_TIMING_GATE_2CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_2CFDFAIL_MAX))){
                filehits++;
                //if (filehits % 1000 == 0)
                //  cout << "Coincident hit " << filehits << " of " << entriesPerFile << ", " << 100 * filehits / entriesPerFile << "% complete" << endl;
              }
            }
          }
        }
        
        if((hit1E >= eLow)&&(hit1E <= eHigh)){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            if(noABHitInd2 != noABHitInd){
              if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                continue; //skip pileup hit
              }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                continue; //skip non-pileup hit
              }

              const double hit2E = offset + sortedEvt.noABHit[noABHitInd2].energy*gain + sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy*quad;

              //if(((sortedEvt.noABHit[noABHitInd].core & 63U)/4)!=((sortedEvt.noABHit[noABHitInd2].core & 63U)/4)){ //try to reduce crosstalk... doesn't seem to do anything regarding sum peak shapes, but seems to align time-random with singles data
              Double_t tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
              uint8_t numCFDFail = 0;
              if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(((numCFDFail == 0)&&(tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX))||((numCFDFail == 1)&&(tDiff >= COINC_TIMING_GATE_1CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_1CFDFAIL_MAX))||((numCFDFail == 2)&&(tDiff >= COINC_TIMING_GATE_2CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_2CFDFAIL_MAX))){
                int eGamma2 = (int)(hit2E/keVPerBin);
                if(eGamma2>=0 && eGamma2<S32K){
                  mcaOut[SP_GATED][eGamma2]++; //fill true coincidence histogram
                }
              }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                //time random
                int eGamma2 = (int)(hit2E/keVPerBin);
                //double eGamma2 = ( (timeRandomOffsetFactor*hit2E*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                if(((int)eGamma2)>=0 && ((int)eGamma2)<S32K){
                  mcaOut[SP_TR_GATED][((int)eGamma2)]++; //fill time-random 'coincidence' histogram
                }
              }
            }
            //}
          }
          //shouldn't break, what if there are 2 gammas in the gate?
          //break;
          
        }
      }

      //make gated summing histograms
      //here we need a triple coincidence, so for each unique set of 3 gammas (A,B,C)
      //we need to determine if any fall within the energy gate, and if they do, whether the other
      //2 gammas are at 180 degrees
      for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
          continue; //skip pileup hit
        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
          continue; //skip non-pileup hit
        }
        const double hit1E = offset + sortedEvt.noABHit[noABHitInd].energy*gain + sortedEvt.noABHit[noABHitInd].energy*sortedEvt.noABHit[noABHitInd].energy*quad;

        if((hit1E >= eLow)&&(hit1E <= eHigh)){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

            if(noABHitInd2 != noABHitInd){
              if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                continue; //skip pileup hit
              }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                continue; //skip non-pileup hit
              }

              uint8_t numCFDFail = 0;
              if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }
              if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                numCFDFail++;
              }

              const double hit2E = offset + sortedEvt.noABHit[noABHitInd2].energy*gain + sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy*quad;

              for(int noABHitInd3 = noABHitInd2+1; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                if(noABHitInd3 == noABHitInd){
                  continue;
                }
                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
                  continue; //skip pileup hit
                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7)))){
                  continue; //skip non-pileup hit
                }

                if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                  //2nd hit and 3rd hit are a unique pair that are 180 degrees apart
                  //now we need to know which of these hits comes first in time

                  int firstSumHitInd = noABHitInd2;
                  int secondSumHitInd = noABHitInd3;
                  if(sortedEvt.noABHit[noABHitInd3].tsDiff < sortedEvt.noABHit[noABHitInd2].tsDiff){
                    firstSumHitInd = noABHitInd3;
                    secondSumHitInd = noABHitInd2;
                  }

                  //Now check the sum timing condition
                  Double_t tDiffSum = (sortedEvt.noABHit[secondSumHitInd].tsDiff - sortedEvt.noABHit[firstSumHitInd].tsDiff)*10.0;
                  if((tDiffSum >= SUM_TIMING_GATE_MIN)&&(tDiffSum <= SUM_TIMING_GATE_MAX)){ //timing condition (sum)

                    const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;

                    //In a true sum event in the DAQ, only the time of the first sum hit matters for evaluating coincidences,
                    //because the 2nd sum hit falls within the programmable deadtime of the first sum
                    //hit and just adds to the energy of the first sum hit. So all coincidences with a gated gamma 
                    //should be evaluated using the time of the first sum hit only

                    const Double_t tDiff = fabs(noABHitTime(&sortedEvt,firstSumHitInd) - noABHitTime(&sortedEvt,noABHitInd));
                    
                    if(((numCFDFail == 0)&&(tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX))||((numCFDFail == 1)&&(tDiff >= COINC_TIMING_GATE_1CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_1CFDFAIL_MAX))||((numCFDFail == 2)&&(tDiff >= COINC_TIMING_GATE_2CFDFAIL_MIN)&&(tDiff <= COINC_TIMING_GATE_2CFDFAIL_MAX))){
                      //time coincident summing
                      const int eSumGamma1 = (int)(hit2E/keVPerBin);
                      const int eSumGamma2 = (int)(hit3E/keVPerBin);
                      const int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy)+(sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
                      //int eGammaSum = (int)(correctSumE(hit2E,hit3E,tDiffSum)/keVPerBin);
                      if(eGammaSum>=0 && eGammaSum<S32K){
                        mcaOut[SP_SUMIN][eGammaSum]++; //fill 180 degree sum histogram
                      }
                      if(eSumGamma2>=0 && eSumGamma2<S32K){
                        mcaOut[SP_SUMOUT][eSumGamma2]++; //fill 180 degree projection histogram
                      }
                      if(eSumGamma1>=0 && eSumGamma1<S32K){
                        mcaOut[SP_SUMOUT][eSumGamma1]++; //fill 180 degree projection histogram
                      }
                      break;
                    }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                      //time random summing
                      const int eSumGamma1 = (int)(hit2E/keVPerBin);
                      const int eSumGamma2 = (int)(hit3E/keVPerBin);
                      const int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy)+(sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
                      //int eGammaSum = (int)(correctSumE(hit2E,hit3E,tDiffSum)/keVPerBin);
                      if(eGammaSum>=0 && eGammaSum<S32K){
                        mcaOut[SP_TR_SUMIN][eGammaSum]++; //fill 180 degree sum histogram
                      }
                      if(eSumGamma2>=0 && eSumGamma2<S32K){
                        mcaOut[SP_TR_SUMOUT][eSumGamma2]++; //fill 180 degree projection histogram
                      }
                      if(eSumGamma1>=0 && eSumGamma1<S32K){
                        mcaOut[SP_TR_SUMOUT][eSumGamma1]++; //fill 180 degree projection histogram
                      }
                      break;
                    }

                  }

                }

              }
              
            }
          }
        }
      }
      

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
  double offset = 0.0;
  double gain = 1.0;
  double quad = 0.0;
  printf("Starting EEGamma_noAB_mca_SMOL_splitruns\n");

  if((argc < 6)||(argc > 11)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_noAB_mca_SMOL_splitruns smol_file_list EGateLow EGateHigh output_dmca_file_prefix events_per_split keV_per_bin discard_pileup offset gain quad" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
    cout << "  *offset*, *gain*, and *quad* are parameters to (re)calibrate the SMOL tree data by. If not specified, these will default to values of 0, 1, and 0 (ie. no change in calibration)." << endl;
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
        discardPileup = atoi(argv[7]);
        if(discardPileup > 2){
          cout << "ERROR: Invalid value for discard_pileup!" << endl;
          cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
          return 0;
        }
        if(argc >= 11){
          offset = atof(argv[8]);
          gain = atof(argv[9]);
          quad = atof(argv[10]);
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
    }else if(discardPileup == 2){
      printf("Only taking pileup hits.\n");
    }
    if(argc == 11){
      printf("Recalibrating with offset = %f, gain = %f, quad = %.15f\n",offset,gain,quad);
    }

    //construct 180 degree summing hit map
    memset(hitMap180deg,0,sizeof(hitMap180deg));
    for(uint8_t i=0;i<64;i++){ //first core
      for(uint8_t j=0;j<64;j++){ //coinc core
        if(i!=j){
          if(getGeVector(i,0,1).Angle(getGeVector(j,0,1))*180.0/PI > 175.0){ //same effect for any value down to 165 degrees
            hitMap180deg[i][j] = 1;
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
          numSepEvts += mysort->SortData(str, outfile, eLow, eHigh, keVPerBin, discardPileup, entriesPerFile, offset, gain, quad);
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
