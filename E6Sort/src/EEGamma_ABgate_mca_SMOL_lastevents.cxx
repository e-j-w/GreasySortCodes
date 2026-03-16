//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_ABgate_mca_SMOL_lastevents_cxx
#include "common.cxx"
#include "EEGamma_ABgate_mca_SMOL_lastevents.h"

using namespace std;

void WriteData(const char* outName, const uint8_t subsetNum, const uint8_t gateNum){

  if(gateNum < MAX_NUM_GATES){
    //duplicate singles data
    if(gateNum > 0){
      memcpy(&mcaOut[gateNum][subsetNum][SP_SINGLES],&mcaOut[0][subsetNum][SP_SINGLES],sizeof(mcaOut[0][subsetNum][SP_SINGLES]));
      memcpy(&mcaOut[gateNum][subsetNum][SP_SINGLES_SUMIN],&mcaOut[0][subsetNum][SP_SINGLES_SUMIN],sizeof(mcaOut[0][subsetNum][SP_SINGLES_SUMIN]));
      memcpy(&mcaOut[gateNum][subsetNum][SP_SINGLES_SUMOUT],&mcaOut[0][subsetNum][SP_SINGLES_SUMOUT],sizeof(mcaOut[0][subsetNum][SP_SINGLES_SUMOUT]));
    }

    printf("Writing subset %u gate %u histogram to: %s\n",subsetNum,gateNum,outName);

    FILE *out;
    if((out = fopen(outName, "w")) == NULL){ //open the file
      printf("ERROR: Cannot open the output file: %s\n",outName);
      return;
    }else{
      fwrite(&mcaOut[gateNum][subsetNum],sizeof(mcaOut[gateNum][subsetNum]),1,out);
      fclose(out);
    }
  }else if(gateNum == 255){
    //special case where there are no gates
    //only write out the singles data

    printf("Writing subset %u singles histogram to: %s\n",subsetNum,outName);

    FILE *out;
    if((out = fopen(outName, "w")) == NULL){ //open the file
      printf("ERROR: Cannot open the output file: %s\n",outName);
      return;
    }else{
      fwrite(&mcaOut[0][subsetNum][SP_SINGLES],sizeof(mcaOut[0][subsetNum][SP_SINGLES]),1,out);
      fwrite(&mcaOut[0][subsetNum][SP_SINGLES_SUMOUT],sizeof(mcaOut[0][subsetNum][SP_SINGLES_SUMOUT]),1,out);
      fwrite(&mcaOut[0][subsetNum][SP_SINGLES_SUMIN],sizeof(mcaOut[0][subsetNum][SP_SINGLES_SUMIN]),1,out);
      fclose(out);
    }
  }
  

}

uint8_t evalCoincGateCFD(const uint8_t numCFDFail, const double tDiff){
  if((numCFDFail == 0)&&(tDiff >= coincGateMin)&&(tDiff <= coincGateMax)){
    return 1; //passes gate
  }else if((numCFDFail == 1)&&(tDiff >= coincGate1CFDFailMin)&&(tDiff <= coincGate1CFDFailMax)){
    return 1; //passes gate
  }else if((numCFDFail == 2)&&(tDiff >= coincGate2CFDFailMin)&&(tDiff <= coincGate2CFDFailMax)){
    return 1; //passes gate
  }else{
    return 0;
  }
}

uint8_t evalSumCorrGateCFD(const uint8_t numCFDFail, const double tDiff){
  if((numCFDFail == 0)&&(tDiff >= sumGateCFDMin)&&(tDiff <= sumGateCFDMax)){
    return 1; //passes gate
  }else if((numCFDFail == 1)&&(tDiff >= sumGate1CFDFailMin)&&(tDiff <= sumGate1CFDFailMax)){
    return 1; //passes gate
  }else if((numCFDFail == 2)&&(tDiff >= sumGate2CFDFailMin)&&(tDiff <= sumGate2CFDFailMax)){
    return 1; //passes gate
  }else{
    return 0;
  }
}

void fillSp(const uint8_t enGate, const uint8_t spInd, const int spBin){
  if((spBin >= 0)&&(spBin < S32K)){
    for(uint8_t i=0; i<numPctToSort; i++){
      if(sortingSubset & (1U << i)){
        mcaOut[enGate][i][spInd][spBin]++;
      }
    }
  }
}

uint64_t getNumEntriesInFile(const char *sfile){
  FILE *inp = fopen(sfile, "rb");
  if(inp == NULL){
    printf("ERROR: couldn't open file %s\n",sfile);
  }
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  fclose(inp);
  return sentries;
}

void SortData(const char *sfile,  
                  const double keVPerBin, const uint8_t discardPileup, 
                  const double offset, const double gain, 
                  const double quad){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  uint64_t pileupCtrs[16];
  fread(&sentries,sizeof(uint64_t),1,inp);
  uint64_t smolVersion = (uint64_t)(sentries >> 48);
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  if((totalEntriesRead + sentries) < (totalEntriesInFileList - maxEvtsToSort)){
    //nothing to sort in this file, move on to the next one
    printf("Skipping file since no events will be sorted.\n");
    totalEntriesRead += sentries;
    return;
  }

  if(smolVersion > 0){
    fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
    //printf("\nNumber of hits of each pileup type:\n");
    uint64_t totalHits = 0;
    for(uint8_t i=0; i<16; i++){
      //printf("Pileup type %2u: %lu\n",i,pileupCtrs[i]);
      totalHits += pileupCtrs[i];
    }
    //printf("Total hits:     %lu\n",totalHits);
    long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
    printf("Fraction of hits with pileup type 1 (no pileup): %Lf\n",frac);
  }

  uint64_t startEntry = 0;
  if(totalEntriesRead < (totalEntriesInFileList - maxEvtsToSort)){
    startEntry = (totalEntriesInFileList - maxEvtsToSort) - totalEntriesRead;
  }

  sorted_evt sortedEvt;

  //allow decimation of sorted events (for debugging/tuning)
  Long64_t increment = 1;
  if(increment == 1){
    printf("\nSorting events (skipping %lu)...\n",startEntry);
  }else if(increment > 0){
    printf("\nSorting every %i events (skipping %lu)...\n",increment,startEntry);
  }else{
    increment = 1;
  }

  for(Long64_t jentry = 0; jentry < sentries; jentry+=increment){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(uint8_t i=0; i<numPctToSort; i++){
      if(!(sortingSubset & (1U << i))){
        if((totalEntriesInFileList - totalEntriesRead) <= evtsToSort[i]){
          sortingSubset |= (1U << i); //flag the subset of data to be sorted
        }
      }
    }
    totalEntriesRead++;
    

    if(jentry < startEntry){
      continue; //don't sort event
    }

    //construct addback energies and times
    memset(addbackE,0,sizeof(addbackE));
    memset(maxABHitE,0,sizeof(maxABHitE));
    memset(addbackNumCFDFail,0,sizeof(addbackNumCFDFail));
    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){
      addbackT[ABpos] = -1.0; //default value
      addbackTS[ABpos] = 255U; //default value
    }
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      
      if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
        continue; //skip pileup hit
      }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
        continue; //skip non-pileup hit
      }

      int ABHitPos = (sortedEvt.noABHit[noABHitInd].core & 63U)/4;

      if(ABHitPos < NGRIFPOS){

        double ABhitE = offset + sortedEvt.noABHit[noABHitInd].energy*gain + sortedEvt.noABHit[noABHitInd].energy*sortedEvt.noABHit[noABHitInd].energy*quad;

        //check timing criteria
        if(addbackT[ABHitPos] >= 0.0){
          //there are hits in this clover
          if(fabs(((double)sortedEvt.noABHit[noABHitInd].tsDiff) - ((double)addbackTS[ABHitPos]))*10.0 > ADDBACK_TIMING_GATE){
            //hit not in time coincidence with other hits
            if(ABhitE > maxABHitE[ABHitPos]){
              //higher energy, not time coincident
              //make this hit the new hit
              addbackT[ABHitPos] = noABHitTime(&sortedEvt,noABHitInd);
              addbackTS[ABHitPos] = sortedEvt.noABHit[noABHitInd].tsDiff;
              if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
                addbackNumCFDFail[ABHitPos]=1;
              }
              maxABHitE[ABHitPos] = ABhitE;
              addbackE[ABHitPos] = ABhitE;
              //got to the next hit
              continue;
            }else{
              //lower energy, not time coincident
              //skip this hit
              continue;
            }
          }
        }
        
        //only get here if there are no hits in the clover, or if there
        //is a time coincident hit
        if(ABhitE > maxABHitE[ABHitPos]){
          addbackT[ABHitPos] = noABHitTime(&sortedEvt,noABHitInd);
          addbackTS[ABHitPos] = sortedEvt.noABHit[noABHitInd].tsDiff;
          if(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)){
            addbackNumCFDFail[ABHitPos]=1;
          }
          maxABHitE[ABHitPos] = ABhitE;
        }
        addbackE[ABHitPos] += ABhitE; 
      }
    }


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
      fillSp(0,SP_SINGLES,singlesE); //fill only for the first energy gate, singles spectra will be duplicated across all energy gates later
      //make summing histograms (for singles)
      for(int noABHitInd3 = noABHitInd+1; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
          continue; //skip pileup hit
        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7)))){
          continue; //skip non-pileup hit
        }
        if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd].core & 63U] != 0){
          Double_t tDiffSum = fabs(sortedEvt.noABHit[noABHitInd3].tsDiff - sortedEvt.noABHit[noABHitInd].tsDiff);
          //printf("tDiffSum: %f\n",tDiffSum);
          if((tDiffSum >= sumGateMin)&&(tDiffSum <= sumGateMax)){ //timing condition (sum)
            //printf("Passed timing condition.\n");

            const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;

            int eGamma3 = (int)(hit3E/keVPerBin);
            int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd].energy+sortedEvt.noABHit[noABHitInd3].energy)*(sortedEvt.noABHit[noABHitInd].energy+sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
            //int eGammaSum = (int)(correctSumE(hit1E,hit3E,tDiffSum)/keVPerBin);
            
            fillSp(0,SP_SINGLES_SUMIN,eGammaSum); //fill 180 degree sum histogram
            fillSp(0,SP_SINGLES_SUMOUT,eGamma3); //fill 180 degree projection histogram
            fillSp(0,SP_SINGLES_SUMOUT,singlesE); //fill 180 degree projection histogram
            
          }
        }
      }
    }

    //look for coincidences with addback hits
    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){

      if(addbackT[ABpos] < 0.0){
        //no addback hit in this clover
        continue;
      }
      
      for(uint8_t i=0; i<numEGates; i++){
        if((addbackE[ABpos] >= gateELow[i])&&(addbackE[ABpos] <= gateEHigh[i])){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            const int ABpos2 = (sortedEvt.noABHit[noABHitInd2].core & 63U)/4;
            if(ABpos2 != ABpos){

              if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                continue; //skip pileup hit
              }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                continue; //skip non-pileup hit
              }

              const double hit2E = offset + sortedEvt.noABHit[noABHitInd2].energy*gain + sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy*quad;

              //if(((sortedEvt.noABHit[noABHitInd].core & 63U)/4)!=((sortedEvt.noABHit[noABHitInd2].core & 63U)/4)){ //try to reduce crosstalk... doesn't seem to do anything regarding sum peak shapes, but seems to align time-random with singles data
              {
                const Double_t tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd2) - addbackT[ABpos]);
                uint8_t numCFDFail = 0;
                if(addbackNumCFDFail[ABpos]){
                  numCFDFail++;
                }
                if(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6)){
                  numCFDFail++;
                }
                if(evalCoincGateCFD(numCFDFail,tDiff)){
                  const int eGamma2 = (int)(hit2E/keVPerBin);
                  fillSp(i,SP_GATED,eGamma2); //fill true coincidence histogram
                }else if((tDiff >= tRandGateMin)&&(tDiff <= tRandGateMax)){
                  //time random
                  const int eGamma2 = (int)(hit2E/keVPerBin);
                  //double eGamma2 = ( (timeRandomOffsetFactor*hit2E*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                  fillSp(i,SP_TR_GATED,((int)eGamma2)); //fill time-random 'coincidence' histogram
                }
              }
              {
                const Double_t tDiffTS = fabs(sortedEvt.noABHit[noABHitInd2].tsDiff - addbackTS[ABpos]);
                if((tDiffTS >= leCoincGateMin)&&(tDiffTS <= leCoincGateMax)){
                  const int eGamma2 = (int)(hit2E/keVPerBin);
                  fillSp(i,SP_LE_GATED,eGamma2); //fill true coincidence histogram
                }else if((tDiffTS >= leTRandGateMin)&&(tDiffTS <= leTRandGateMax)){
                  //time random
                  const int eGamma2 = (int)(hit2E/keVPerBin);
                  //double eGamma2 = ( (timeRandomOffsetFactor*hit2E*(1.0+rand_sym_dbl(timeRandomWidthInfl)))/keVPerBin );
                  fillSp(i,SP_LE_TR_GATED,((int)eGamma2)); //fill time-random 'coincidence' histogram
                }
              }

              //make gated summing histograms
              //here we need a triple coincidence, so for each unique set of 3 gammas (A,B,C)
              //we need to determine if any fall within the energy gate, and if they do, whether the other
              //2 gammas are at 180 degrees
              for(int noABHitInd3 = noABHitInd2+1; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){

                //getting rid of both these conditions over-estimates summing,
                //but getting rid of neither under-estimates summing
                const int ABpos3 = (sortedEvt.noABHit[noABHitInd3].core & 63U)/4;
                if(ABpos3 == ABpos){
                  //const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;
                  //if((hit3E >= gateELow[i])&&(hit3E <= gateEHigh[i])){
                    continue;
                  //}
                }
                
                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
                  continue; //skip pileup hit
                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7)))){
                  continue; //skip non-pileup hit
                }

                if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                  //2nd hit and 3rd hit are a unique pair that are 180 degrees apart
                  //now we need to know which of these hits comes first in time,
                  //since that gives the DAQ time of the equivalent 0 degree sum hit

                  int firstSumHitInd = noABHitInd2;
                  int secondSumHitInd = noABHitInd3;
                  if(sortedEvt.noABHit[noABHitInd3].tsDiff < sortedEvt.noABHit[noABHitInd2].tsDiff){
                    firstSumHitInd = noABHitInd3;
                    secondSumHitInd = noABHitInd2;
                  }
                  int firstSumHitIndCFD = noABHitInd2;
                  int secondSumHitIndCFD = noABHitInd3;
                  if(noABHitTime(&sortedEvt,noABHitInd3) < noABHitTime(&sortedEvt,noABHitInd2)){
                    firstSumHitIndCFD = noABHitInd3;
                    secondSumHitIndCFD = noABHitInd2;
                  }

                  //Now check the sum timing condition
                  Double_t tDiffSum = sortedEvt.noABHit[secondSumHitInd].tsDiff - sortedEvt.noABHit[firstSumHitInd].tsDiff;
                  //printf("tDiffSum: %f\n",tDiffSum);
                  if((tDiffSum >= sumGateMin)&&(tDiffSum <= sumGateMax)){ //timing condition (sum)

                    const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;
                    const int eSumGamma1 = (int)(hit2E/keVPerBin);
                    const int eSumGamma2 = (int)(hit3E/keVPerBin);
                    const int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*(sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
                    //int eGammaSum = (int)(correctSumE(hit2E,hit3E,tDiffSum)/keVPerBin);

                    //In a true sum event in the DAQ, only the time of the first sum hit matters for evaluating coincidences,
                    //because the 2nd sum hit falls within the programmable deadtime of the first sum
                    //hit and just adds to the energy of the first sum hit. So all coincidences with a gated gamma 
                    //should be evaluated using the time of the first sum hit only

                    {

                      const Double_t tDiff = fabs(noABHitTime(&sortedEvt,firstSumHitInd) - addbackT[ABpos]);
                      uint8_t numCFDFail = 0;
                      if(addbackNumCFDFail[ABpos]){
                        numCFDFail++;
                      }
                      if(sortedEvt.noABHit[firstSumHitInd].core & ((uint8_t)1 << 6)){
                        numCFDFail++;
                      }

                      if(evalCoincGateCFD(numCFDFail,tDiff)){
                        //time coincident summing
                        fillSp(i,SP_SUMIN,eGammaSum); //fill 180 degree sum histogram
                        fillSp(i,SP_SUMOUT,eSumGamma2); //fill 180 degree projection histogram
                        fillSp(i,SP_SUMOUT,eSumGamma1); //fill 180 degree projection histogram
                      }else if((tDiff >= tRandGateMin)&&(tDiff <= tRandGateMax)){
                        //time random summing
                        fillSp(i,SP_TR_SUMIN,eGammaSum); //fill 180 degree sum histogram
                        fillSp(i,SP_TR_SUMOUT,eSumGamma2); //fill 180 degree projection histogram
                        fillSp(i,SP_TR_SUMOUT,eSumGamma1); //fill 180 degree projection histogram
                      }
                    }
                    
                    {
                      const Double_t tDiffTS = fabs(sortedEvt.noABHit[firstSumHitInd].tsDiff - addbackTS[ABpos]);
                      //printf("tDiffTS: %f (%u %u)\n",tDiffTS,sortedEvt.noABHit[firstSumHitInd].tsDiff,addbackTS[ABpos]);
                      if((tDiffTS >= leCoincGateMin)&&(tDiffTS <= leCoincGateMax)){
                        //time coincident summing
                        fillSp(i,SP_LE_SUMIN,eGammaSum); //fill 180 degree sum histogram
                        fillSp(i,SP_LE_SUMOUT,eSumGamma2); //fill 180 degree projection histogram
                        fillSp(i,SP_LE_SUMOUT,eSumGamma1); //fill 180 degree projection histogram
                      }else if((tDiffTS >= leTRandGateMin)&&(tDiffTS <= leTRandGateMax)){
                        //time random summing
                        fillSp(i,SP_LE_TR_SUMIN,eGammaSum); //fill 180 degree sum histogram
                        fillSp(i,SP_LE_TR_SUMOUT,eSumGamma2); //fill 180 degree projection histogram
                        fillSp(i,SP_LE_TR_SUMOUT,eSumGamma1); //fill 180 degree projection histogram
                      }
                    }

                  }

                  //check the sum timing condition, using CFD timing
                  Double_t tDiffSumCFD = noABHitTime(&sortedEvt,secondSumHitIndCFD) - noABHitTime(&sortedEvt,firstSumHitIndCFD);
                  uint8_t numCFDFailSum = 0;
                  if(sortedEvt.noABHit[firstSumHitIndCFD].core & ((uint8_t)1 << 6)){
                    numCFDFailSum++;
                  }
                  if(sortedEvt.noABHit[secondSumHitIndCFD].core & ((uint8_t)1 << 6)){
                    numCFDFailSum++;
                  }
                  if(evalSumCorrGateCFD(numCFDFailSum, tDiffSumCFD)){

                    //In a true sum event in the DAQ, only the time of the first sum hit matters for evaluating coincidences,
                    //because the 2nd sum hit falls within the programmable deadtime of the first sum
                    //hit and just adds to the energy of the first sum hit. So all coincidences with a gated gamma 
                    //should be evaluated using the time of the first sum hit only
                    {
                      const Double_t tDiff = fabs(noABHitTime(&sortedEvt,firstSumHitIndCFD) - addbackT[ABpos]);
                      uint8_t numCFDFail = 0;
                      if(addbackNumCFDFail[ABpos]){
                        numCFDFail++;
                      }
                      if(sortedEvt.noABHit[firstSumHitIndCFD].core & ((uint8_t)1 << 6)){
                        numCFDFail++;
                      }
                      if(evalCoincGateCFD(numCFDFail,tDiff)){
                        //time coincident summing
                        const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;
                        const int eSumGamma1 = (int)(hit2E/keVPerBin);
                        const int eSumGamma2 = (int)(hit3E/keVPerBin);
                        const int eGammaSum = (int)((offset + (sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*gain + ((sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy)*(sortedEvt.noABHit[noABHitInd2].energy+sortedEvt.noABHit[noABHitInd3].energy))*quad)/keVPerBin);
                        //int eGammaSum = (int)(correctSumE(hit2E,hit3E,tDiffSum)/keVPerBin);
                        fillSp(i,SP_SUMIN_CFD,eGammaSum); //fill 180 degree sum histogram
                        fillSp(i,SP_SUMOUT_CFD,eSumGamma2); //fill 180 degree projection histogram
                        fillSp(i,SP_SUMOUT_CFD,eSumGamma1); //fill 180 degree projection histogram
                      }
                    }

                  }

                }

              }


            }
            //}
          }
          //shouldn't break, what if there are 2 gammas in the gate?
          //break;
          
        }
      }
      
    }

    if (jentry % 90713 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << (jentry-startEntry) << " of " << (sentries-startEntry) << ", " << 100 * (jentry-startEntry) / (sentries-startEntry) << "% complete" << "\r" << flush;

  } // analysis tree

  cout << "Entry " << (sentries-startEntry) << " of " << (sentries-startEntry) << ", 100% complete" << endl;
  
  fclose(inp);
  
}

int main(int argc, char **argv){

  const char *sfile;
  const char *outfile;
  uint8_t discardPileup = 0;
  double keVPerBin = 1.0;
  double lastEvtsPercent = 0.0;
  memset(evtsToSort,0,sizeof(evtsToSort));
  maxEvtsToSort = 0;
  double offset = 0.0;
  double gain = 1.0;
  double quad = 0.0;
  sortingSubset = 0;
  coincGateMin = COINC_TIMING_GATE_MIN;
  coincGateMax = COINC_TIMING_GATE_MAX;
  coincGate1CFDFailMin = COINC_TIMING_GATE_1CFDFAIL_MIN;
  coincGate1CFDFailMax = COINC_TIMING_GATE_1CFDFAIL_MAX;
  coincGate2CFDFailMin = COINC_TIMING_GATE_2CFDFAIL_MIN;
  coincGate2CFDFailMax = COINC_TIMING_GATE_2CFDFAIL_MAX;
  sumGateMin = SUM_TIMING_GATE_MIN;
  sumGateMax = SUM_TIMING_GATE_MAX;
  tRandGateMin = TRANDOM_GATE_MIN;
  tRandGateMax = TRANDOM_GATE_MAX;
  leCoincGateMin = LE_COINC_TIMING_GATE_MIN;
  leCoincGateMax = LE_COINC_TIMING_GATE_MAX;
  leTRandGateMin = LE_TRANDOM_GATE_MIN;
  leTRandGateMax = LE_TRANDOM_GATE_MAX;
  printf("Starting EEGamma_ABgate_mca_SMOL_lastevents\n");

  if(argc < 6){
    cout << "Generates GRIFIFN gated spectra." << endl;
    cout << "Arguments: EEGamma_ABgate_mca_SMOL_lastevents smol_file_list num_gates EGateLow1 EGateHigh1 (EGateLow2 EGateHigh2...) num_sorts percent_of_events_1 (percent_of_events_2 ...) output_dmca_file_prefix keV_per_bin discard_pileup offset gain quad" << endl;
    cout << "  *smol_file* must be a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *percent_of_events_X* specifies the percentage of events at the end of the file list to sort. The intention when writing this was to sort only events at the end of a decay curve." << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
    cout << "  *offset*, *gain*, and *quad* are parameters to (re)calibrate the SMOL tree data by. If not specified, these will default to values of 0, 1, and 0 (ie. no change in calibration)." << endl;
    cout << "  Additional parameters corresponding to timing gates can optionally be specified afterward, in the format:" << endl;
    cout << "    coincGateMin coincGateMax coincGate1CFDFailMin coincGate1CFDFailMax coincGate2CFDFailMin coincGate2CFDFailMax sumGateMin sumGateMax tRandGateMin tRandGateMax leCoincGateMin leCoincGateMax leTRandGateMin leTRandGateMax" << endl;
    return 0;
  }else{
    sfile = argv[1];
    numEGates = atoi(argv[2]);
    uint8_t currentArg = 3;
    if((numEGates >= 0)&&(numEGates <= MAX_NUM_GATES)){
      //valid number of gates
      //check that there are enough arguments
      if(argc < (5 + 2*numEGates)){
        printf("ERROR: not enough arguments for the number of energy gates specified (need %u).\n",(5 + 2*numEGates));
        return 0;
      }
      for(uint8_t i=0; i<numEGates; i++){
        gateELow[i] = atof(argv[currentArg++]);
        gateEHigh[i] = atof(argv[currentArg++]);
      }
    }else{
      printf("ERROR: invalid number of energy gates (%u).\n",numEGates);
      return 0;
    }
    numPctToSort = (uint8_t)(atoi(argv[currentArg++]));
    //printf("Number of subsets to sort: %u.\n",numPctToSort);
    if((numPctToSort > 0)&&(numPctToSort <= MAX_NUM_PCTTOSORT)){
      //valid number of subsets of data to sort
      if(argc < (5 + 2*numEGates + numPctToSort)){
        printf("ERROR: not enough arguments for the number of sorts specified (need %u).\n",(5 + 2*numEGates + numPctToSort));
        return 0;
      }
      for(uint8_t i=0; i<numPctToSort; i++){
        pctToSort[i] = atof(argv[currentArg++]);
        //printf("%f\n",pctToSort[i]);
      }
    }else{
      printf("ERROR: invalid number of sorts specified (%u).\n",numPctToSort);
      return 0;
    }
    outfile = argv[currentArg++];
    //printf("Output filepath: %s.\n",argv[currentArg-1]);
    if(argc > (5 + 2*numEGates + numPctToSort)){
      keVPerBin = atof(argv[currentArg++]);
      if(argc > (6 + 2*numEGates + numPctToSort)){
        discardPileup = atoi(argv[currentArg++]);
        if(discardPileup > 2){
          printf("ERROR: Invalid value for discard_pileup (%s)!\n",argv[currentArg-1]);
          printf("  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits).\n");
          return 0;
        }
        if(argc >= (10 + 2*numEGates)){
          offset = atof(argv[currentArg++]);
          gain = atof(argv[currentArg++]);
          quad = atof(argv[currentArg++]);
        }
        if(argc >= (24 + 2*numEGates)){
          //manually specified timing gates
          coincGateMin = atof(argv[currentArg++]);
          coincGateMax = atof(argv[currentArg++]);
          coincGate1CFDFailMin = atof(argv[currentArg++]);
          coincGate1CFDFailMax = atof(argv[currentArg++]);
          coincGate2CFDFailMin = atof(argv[currentArg++]);
          coincGate2CFDFailMax = atof(argv[currentArg++]);
          sumGateMin = atof(argv[currentArg++]);
          sumGateMax = atof(argv[currentArg++]);
          tRandGateMin = atof(argv[currentArg++]);
          tRandGateMax = atof(argv[currentArg++]);
          leCoincGateMin = atof(argv[currentArg++]);
          leCoincGateMax = atof(argv[currentArg++]);
          leTRandGateMin = atof(argv[currentArg++]);
          leTRandGateMax = atof(argv[currentArg++]);
        }
      }

    }
  }

  for(uint8_t i=0; i<numPctToSort; i++){
    if((pctToSort[i] <= 0.0)||(pctToSort[i] > 100.0)){
      printf("Invalid event percentage to sort (%f)!\nThe event percentage must be greater than > 0 and <= 100.\n",pctToSort[i]);
      return 0;
    }
  }

  if(keVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
    return 0;
  }
  
  for(uint8_t i=0; i<numEGates; i++){
    if(gateELow[i] > gateEHigh[i]){
      //swap energy gate bounds
      double tmp = gateEHigh[i];
      gateEHigh[i] = gateELow[i];
      gateELow[i] = tmp;
    }
  }

  //setup sum correction CFD gates
  //basically, make sure that the gate start and gate width are offset the same as the coincidence gates
  sumGateCFDMin = sumGateMin*10.0;
  sumGateCFDMax = sumGateMax*10.0;
  sumGate1CFDFailMin = (coincGate1CFDFailMin-coincGateMin)+sumGateCFDMin;
  sumGate1CFDFailMax = (fabs(coincGate1CFDFailMax-coincGate1CFDFailMin) - fabs(coincGateMax-coincGateMin)) + fabs(sumGateCFDMax - sumGateCFDMin) + sumGate1CFDFailMin;
  sumGate2CFDFailMin = (coincGate2CFDFailMin-coincGateMin)+sumGateCFDMin;
  sumGate2CFDFailMax = (fabs(coincGate2CFDFailMax-coincGate2CFDFailMin) - fabs(coincGateMax-coincGateMin)) + fabs(sumGateCFDMax - sumGateCFDMin) + sumGate2CFDFailMin;

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum
  numFilesWritten = 0;

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree list file name." << endl;
    return 0;
  }

  if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nEnergy gate(s): [", sfile);
    for(uint8_t i=0; i<numEGates; i++){
      if(i==0){
        printf("%0.2f %0.2f",gateELow[i],gateEHigh[i]);
      }else{
        printf("], [%0.2f %0.2f",gateELow[i],gateEHigh[i]);
      }
    }
    printf("] keV\nOutput file prefix: %s\nPercentage of events to sort: [", outfile);
    for(uint8_t i=0; i<numPctToSort; i++){
      if(i==0){
        printf("%0.2f",pctToSort[i]);
      }else{
        printf("], [%0.2f",pctToSort[i]);
      }
    }
    printf("]\n%0.2f keV per bin\n", keVPerBin);
    printf("Coincidence timing gate: [%0.2f %0.2f] ns\n",coincGateMin,coincGateMax);
    printf("  (with 1 CFD fail: [%0.2f %0.2f] ns)\n",coincGate1CFDFailMin,coincGate1CFDFailMax);
    printf("  (with 2 CFD fails: [%0.2f %0.2f] ns)\n",coincGate2CFDFailMin,coincGate2CFDFailMax);
    printf("Sum timing gate: [%0.2f %0.2f] timestamps\n",sumGateMin,sumGateMax);
    printf("  (using CFD timing: [%0.2f %0.2f] ns)\n",sumGateCFDMin,sumGateCFDMax);
    printf("    (with 1 CFD fail: [%0.2f %0.2f] ns)\n",sumGate1CFDFailMin,sumGate1CFDFailMax);
    printf("    (with 2 CFD fails: [%0.2f %0.2f] ns)\n",sumGate2CFDFailMin,sumGate2CFDFailMax);
    printf("Random timing gate: [%0.2f %0.2f] ns\n",tRandGateMin,tRandGateMax);
    printf("Leading-edge coincidence timing gate: [%0.2f %0.2f] timestamps\n",leCoincGateMin,leCoincGateMax);
    printf("Leading-edge random timing gate: [%0.2f %0.2f] timestamps\n",leTRandGateMin,leTRandGateMax);
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
          /*if(i==1){
            printf("Angle between %u and %u: %f\n",i,j,getGRIFFINVector(i,1).Angle(getGRIFFINVector(j,1))*180.0/PI);
          }*/
          //uses vectors for detectors in forward position - should be valid for back position as well
          if(getGRIFFINVector(i,1).Angle(getGRIFFINVector(j,1))*180.0/PI > 170.0){ //same effect for any value down to 165 degrees
            hitMap180deg[i][j] = 1;
            //printf("Pair %u and %u are at 180 degrees.\n",i,j);
          }
        }
      }
    }
    
    FILE *listfile;
    char str[256];

    totalEntriesRead = 0;
    totalEntriesInFileList = 0;

    //count entries in each file
    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      printf("Determining the total number of events in all trees specified in the file list: %s\n",sfile);
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          printf("Checking file: %s\n",str);
          totalEntriesInFileList += getNumEntriesInFile(str);
        }
      }
      for(uint8_t i=0; i<numPctToSort; i++){
        if(pctToSort[i] >= 100.0){
          evtsToSort[i] = totalEntriesInFileList;
        }else{
          evtsToSort[i] = (uint64_t)(totalEntriesInFileList*pctToSort[i]/100.0);
        }
        if(evtsToSort[i] > maxEvtsToSort){
          maxEvtsToSort = evtsToSort[i];
        }
      }
      
      printf("%lu total events found.\nWill sort: [",totalEntriesInFileList);
      for(uint8_t i=0; i<numPctToSort; i++){
        if(i==0){
          printf("%lu",evtsToSort[i]);
        }else{
          printf("], [%lu",evtsToSort[i]);
        }
      }
      printf("] events.\n");
    }


    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          SortData(str, keVPerBin, discardPileup, 
                   offset, gain, quad);
        }
      }
    }
  }else{
    cout << "ERROR: improper file extension for SMOL tree list (should be .list)." << endl;
    return 0;
  }
  
  char fileName[256];
  for(uint8_t i=0; i<numPctToSort; i++){
    if(numEGates > 0){
      for(uint8_t j=0; j<numEGates; j++){
        int meanGateE = (int)((gateEHigh[j] + gateELow[j])/2.0);
        snprintf(fileName,256,"%s_%.1fpct_%ikeVgate.dmca",outfile,pctToSort[i],meanGateE);
        WriteData(fileName,i,j); //write data to disk
      }
    }else{
      snprintf(fileName,256,"%s_%.1fpct_singles.dmca",outfile,pctToSort[i]);
      WriteData(fileName,i,255); //write data to disk (no energy gates)
    }
  }

  printf("Sorted a total of %lu events, keeping the last [",totalEntriesRead);
  for(uint8_t i=0; i<numPctToSort; i++){
    if(i==0){
      printf("%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
    }else{
      printf("], [%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
    }
  }
  printf("]\n");

  return 0;
}
