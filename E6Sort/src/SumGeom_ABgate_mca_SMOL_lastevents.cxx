//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define SumGeom_ABgate_mca_SMOL_lastevents_cxx
#include "common.cxx"
#include "SumGeom_ABgate_mca_SMOL_lastevents.h"

using namespace std;

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
    if(totalHits > 0){
      long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
      printf("Fraction of hits with pileup type 1 (no pileup): %Lf\n",frac);
    }
  }

  uint64_t startEntry = 0;
  if(totalEntriesRead < (totalEntriesInFileList - maxEvtsToSort)){
    startEntry = (totalEntriesInFileList - maxEvtsToSort) - totalEntriesRead;
  }

  sorted_evt sortedEvt;

  //allow decimation of sorted events (for debugging/tuning)
  uint64_t increment = 1;
  if(increment == 1){
    printf("\nSorting events (skipping %lu)...\n",startEntry);
  }else if(increment > 0){
    printf("\nSorting every %i events (skipping %lu)...\n",increment,startEntry);
  }else{
    increment = 1;
  }

  for(uint64_t jentry = 0; jentry < sentries; jentry+=increment){

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

    //'randomly' sample an addback position to compare to the sum hits
    //while taking into account the efficiency at each clover
    //we do this by taking the first addback hit of an event, then using
    //it to compare against hits in the next event

    if(prevEvtABPos == 255U){
      for(uint8_t ABpos = 0; ABpos < NGRIFPOS; ABpos++){

        if(addbackT[ABpos] < 0.0){
          //no addback hit in this clover
          continue;
        }else if(addbackE[ABpos] >= E_THRESHOLD){ //ignore threshold effects for low energy gammas (we assume any energy gates set are above threshold for all detectors)
          //flag this addback position
          prevEvtABPos = ABpos;
          ABHitEvtNum = jentry;
          break;
        }

        if(ABpos == NGRIFPOS-1){
          //no addback hits found
          prevEvtABPos = 255U; //set invalid position
        }

      }
    }

    //look for coincidences with addback hits
    if(jentry != ABHitEvtNum){ //event mixing
      if(prevEvtABPos != 255U){

        if(prevEvtSumHit1Pos == 255U){
          for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
              continue; //skip pileup hit
            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
              continue; //skip non-pileup hit
            }

            const double hit2E = offset + sortedEvt.noABHit[noABHitInd2].energy*gain + sortedEvt.noABHit[noABHitInd2].energy*sortedEvt.noABHit[noABHitInd2].energy*quad;
            if((hit2E/keVPerBin) >= E_THRESHOLD){ //ignore threshold effects for low energy gammas (we assume any sum peaks being analyzed consist of gammas above threshold for all detectors)
              //flag 2nd hit
              prevEvtSumHit1Pos = sortedEvt.noABHit[noABHitInd2].core & 63U;
              sumHit1EvtNum = jentry;
            }
          }
        }

        if(jentry != sumHit1EvtNum){
          if(prevEvtSumHit1Pos != 255U){
            //evaluate 180 degree summing conditions
            //only look at gammas that have real coincidences, since some positions
            //may not contain any gammas in coincidence with other and therefore
            //wont be available for the 180 degree coincidence summing correction,
            //but will still contain real sum peaks
            //(eg. if the corresponding GRIF-16 has a clock de-sync)
            if(sortedEvt.header.numNoABHits > 1){ //enforce coincidence condition
              for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                
                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7))){
                  continue; //skip pileup hit
                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd3].core & ((uint8_t)(1) << 7)))){
                  continue; //skip non-pileup hit
                }

                if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][prevEvtSumHit1Pos] != 0){
                  //2nd hit and 3rd hit are a unique pair that are 180 degrees apart
                  //(assume that the 2nd hit is the first of the pair, and that they are in coincidence)

                  const double hit3E = offset + sortedEvt.noABHit[noABHitInd3].energy*gain + sortedEvt.noABHit[noABHitInd3].energy*sortedEvt.noABHit[noABHitInd3].energy*quad;
                  if((hit3E/keVPerBin) >= E_THRESHOLD){ //ignore threshold effects for low energy gammas (we assume any sum peaks being analyzed consist of gammas above threshold for all detectors)
                    //flag 3rd hit
                    prevEvtSumHit2Pos = sortedEvt.noABHit[noABHitInd3].core & 63U;
                  }
                  
                }
                //}
              }
            }
            
          }
        }

        //if all 3 hits have been flagged, evaluate what type
        //of summing/correction they contribute to 
        if(prevEvtSumHit2Pos != 255U){
          //now check whether the 180 degree coindident hit conflicts with the original addback gate
          //and flag it if so
          const int ABpos2 = (prevEvtSumHit1Pos)/4;
          const int ABpos3 = (prevEvtSumHit2Pos)/4;
          if((ABpos3 != prevEvtABPos)&&(ABpos2 != prevEvtABPos)){
            //in this case both a real sum peak and a 180 correction would be visible
            //since neither of the 180 degree gammas is in the clover where the coincident hit occured
            for(uint8_t i=0; i<numPctToSort; i++){
              if(sortingSubset & (1U << i)){
                numEvtCoinc[i]++; //real sum peak
                numEvtCoincSum[i]++; //180 sum coincidence
              }
            }
          }else if(ABpos2 != prevEvtABPos){
            //in this case a real sum peak would be visible, but not a 180 degree correction
            //since only one of the 180 degree gammas is in the clover where the coincident hit occured
            for(uint8_t i=0; i<numPctToSort; i++){
              if(sortingSubset & (1U << i)){
                numEvtCoinc[i]++; //real sum peak
              }
            }
          }

          //reset flags
          prevEvtABPos = 255U;
          prevEvtSumHit1Pos = 255U;
          prevEvtSumHit2Pos = 255U;
        }else if(jentry > (ABHitEvtNum+EVT_MIX_SEARCH_DEPTH)){
          //reached the end of the event mixing search
          //but the 3rd hit was not seen
          //check whether the 2nd hit contributes to summing

          if(prevEvtSumHit1Pos != 255U){
            const int ABpos2 = (prevEvtSumHit1Pos)/4;
            if(ABpos2 != prevEvtABPos){
              //in this case a real sum peak would be visible, but not a 180 degree correction
              //since the only other gamma isn't in the clover where the coincident hit occured
              //and 180 degree coincident gammas apparently don't exist
              for(uint8_t i=0; i<numPctToSort; i++){
                if(sortingSubset & (1U << i)){
                  numEvtCoinc[i]++; //real sum peak
                }
              }
            }
          }

          //reset flags
          prevEvtABPos = 255U;
          prevEvtSumHit1Pos = 255U;
          prevEvtSumHit2Pos = 255U;
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
  uint8_t forwardPos = 1;
  sortingSubset = 0;
  printf("Starting SumGeom_ABgate_mca_SMOL_lastevents\n");

  if(argc < 5){
    cout << "Determines the geometrical effect on the 180 degree summing correction, when using an addback gate." << endl;
    cout << "Arguments: SumGeom_ABgate_mca_SMOL_lastevents smol_file_list output_file forward_pos num_sorts percent_of_events_1 (percent_of_events_2 ...) keV_per_bin discard_pileup offset gain quad" << endl;
    cout << "  *smol_file* must be a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *output_file* is a text file that the geometric correction will be written to." << endl;
    cout << "  *percent_of_events_X* specifies the percentage of events at the end of the file list to sort. The intention when writing this was to sort only events at the end of a decay curve." << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
    cout << "  *offset*, *gain*, and *quad* are parameters to (re)calibrate the SMOL tree data by. If not specified, these will default to values of 0, 1, and 0 (ie. no change in calibration)." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    forwardPos = (uint8_t)atoi(argv[3]);
    uint8_t currentArg = 4;
    numPctToSort = (uint8_t)(atoi(argv[currentArg++]));
    //printf("Number of subsets to sort: %u.\n",numPctToSort);
    if((numPctToSort > 0)&&(numPctToSort <= MAX_NUM_PCTTOSORT)){
      //valid number of subsets of data to sort
      if(argc < (5 + numPctToSort)){
        printf("ERROR: not enough arguments for the number of sorts specified (need %u).\n",(5 + numPctToSort));
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
    //printf("Output filepath: %s.\n",argv[currentArg-1]);
    if(argc > (5 + numPctToSort)){
      keVPerBin = atof(argv[currentArg++]);
      if(argc > (6 + numPctToSort)){
        discardPileup = atoi(argv[currentArg++]);
        if(discardPileup > 2){
          printf("ERROR: Invalid value for discard_pileup (%s)!\n",argv[currentArg-1]);
          printf("  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits).\n");
          return 0;
        }
        if(argc >= (10 + numPctToSort)){
          offset = atof(argv[currentArg++]);
          gain = atof(argv[currentArg++]);
          quad = atof(argv[currentArg++]);
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

  //initialize counters
  memset(numEvtCoinc,0,sizeof(numEvtCoinc));
  memset(numEvtCoincSum,0,sizeof(numEvtCoincSum));
  prevEvtABPos = 255U;
  prevEvtSumHit1Pos = 255U;
  prevEvtSumHit2Pos = 255U;

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree list file name." << endl;
    return 0;
  }

  if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\n", sfile);
    printf("Output file: %s\n", outfile);
    if(forwardPos == 1){
      printf("GRIFFIN at 110 mm\n");
    }else if(forwardPos == 0){
      printf("GRIFFIN at 145 mm\n");
    }else{
      printf("ERROR: invalid GRIFFIN position!\n");
      return 0;
    }
    printf("Percentage of events to sort: [");
    for(uint8_t i=0; i<numPctToSort; i++){
      if(i==0){
        printf("%0.2f",pctToSort[i]);
      }else{
        printf("], [%0.2f",pctToSort[i]);
      }
    }
    printf("]\n%0.2f keV per bin\n", keVPerBin);
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
            printf("Angle between %u and %u: %f\n",i,j,getGRIFFINVector(i,forwardPos).Angle(getGRIFFINVector(j,forwardPos))*180.0/PI);
          }*/
          if(getGRIFFINVector(i,forwardPos).Angle(getGRIFFINVector(j,forwardPos))*180.0/PI > 170.0){ //same effect for any value down to 165 degrees
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

  printf("\nSorted a total of %lu events, keeping the last [",totalEntriesRead);
  for(uint8_t i=0; i<numPctToSort; i++){
    if(i==0){
      printf("%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
    }else{
      printf("], [%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
    }
  }
  printf("]\n\n");

  for(uint8_t i=0; i<numPctToSort; i++){
    printf("Sort %: %6.2f, Real geometric summing: %10lu, 180 degree summing: %10lu, ratio: %f\n",pctToSort[i],numEvtCoinc[i],numEvtCoincSum[i],((double)(numEvtCoinc[i]))/((double)(numEvtCoincSum[i])));
  }
  printf("For reference:\n");
  printf(" 1 clover missing: %f\n",15.0/14.0);

  FILE *out;
  if((out = fopen(outfile, "w")) == NULL){ //open the file
    printf("ERROR: Cannot open the output file: %s\n",outfile);
    return 0;
  }else{
    fprintf(out,"Sorted a total of %lu events, keeping the last [",totalEntriesRead);
    for(uint8_t i=0; i<numPctToSort; i++){
      if(i==0){
        fprintf(out,"%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
      }else{
        fprintf(out,"], [%lu (%f %%)",evtsToSort[i],100.0*(evtsToSort[i]/(1.0*totalEntriesRead)));
      }
    }
    fprintf(out,"]\n\n");

    for(uint8_t i=0; i<numPctToSort; i++){
      fprintf(out,"Sort %: %6.2f, Real geometric summing: %10lu, 180 degree summing: %10lu, ratio: %f\n",pctToSort[i],numEvtCoinc[i],numEvtCoincSum[i],((double)(numEvtCoinc[i]))/((double)(numEvtCoincSum[i])));
    }
    fprintf(out,"For reference:\n");
    fprintf(out," 1 clover missing: %f\n",15.0/14.0);
    fclose(out);
  }

  return 0;
}
