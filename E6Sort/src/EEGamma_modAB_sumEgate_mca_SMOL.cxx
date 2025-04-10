//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_modAB_sumEgate_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_modAB_sumEgate_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

uint8_t ABHitMapping[64]; //arrays specifying which hits correspond to which addback hits
double addbackE[64];
double addbackT[64];

void EEGamma_modAB_sumEgate_mca_SMOL::WriteData(const char* outName){

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

uint64_t EEGamma_modAB_sumEgate_mca_SMOL::SortData(const char *sfile, const double eLow, const double eHigh, const double ABRad, const double keVPerBin){

    FILE *inp = fopen(sfile, "rb");
    printf("\nFile %s opened\n", sfile);
    
    uint64_t sentries = 0U;
    fread(&sentries,sizeof(uint64_t),1,inp);
    sorted_evt sortedEvt;
    uint8_t footerVal;

    //construct 180 degree summing hit map
    memset(hitMap180deg,0,sizeof(hitMap180deg));
    for(uint8_t i=0;i<64;i++){ //first core
        for(uint8_t j=0;j<64;j++){ //coinc core
            if(i!=j){
                if(getGeVector(i,0,1).Angle(getGeVector(j,0,1))*180.0/PI > 175.0){
                    hitMap180deg[i][j] = 1;
                    continue; //check the next coinc core
                }
                for(uint8_t k=0;k<64;k++){ //other core next to the first core, which could be addback'd with ti
                    if((k!=i)&&(k!=j)){
                        if(getGeHitDistance(i,0,k,0,1) < ABRad){
                            if(getGeVector(k,0,1).Angle(getGeVector(j,0,1))*180.0/PI > 175.0){
                                hitMap180deg[i][j] = 1;
                                break; //check the next coinc core
                            }
                        }
                    }
                }
            }
        }
    }

    printf("Sorting events...\n");
    for(Long64_t jentry = 0; jentry < sentries; jentry++){

        //read event
        if(readSMOLEvent(inp,&sortedEvt)==0){
            cout << "ERROR: bad event data in entry " << jentry << "." << endl;
            exit(-1);
        }

        if(sortedEvt.header.numNoABHits < 2){
            continue;
        }

        //reset flags
        uint64_t ABHitBuildFlags = 0;
        uint64_t specFillFlags[2]; //flags specifying which hits have been used to fill spectra
        memset(specFillFlags,0,sizeof(specFillFlags));
        uint8_t numABHitsBuilt = 0;
        memset(ABHitMapping,0,sizeof(ABHitMapping));
        memset(addbackE,0,sizeof(addbackE));
        memset(addbackT,0,sizeof(addbackT));

        //build initial addback hits, using the gate addback radius,
        //and using ABHitMapping to track which hits are grouped together
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            if(noABHitInd < 64){
                uint8_t abHitBuilt = 0;
                for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                    if(noABHitInd2 < 64){
                        if(noABHitInd2 != noABHitInd){
                            //check if hits are in neighbouring crystals
                            if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < ABRad){ //FORWARD POSITION (11 cm)
                                double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                if(tDiff <= ADDBACK_TIMING_GATE){ //timing condition
                                    if(!(ABHitBuildFlags & (1UL << noABHitInd))){
                                        //first hit not yet flagged
                                        if(!(ABHitBuildFlags & (1UL << noABHitInd2))){
                                            //neither hit flagged
                                            //build initial addback hit
                                            addbackT[numABHitsBuilt] = noABHitTime(&sortedEvt,noABHitInd);
                                            if(!(ABHitBuildFlags & (1UL << noABHitInd))){
                                                addbackE[numABHitsBuilt] += sortedEvt.noABHit[noABHitInd].energy;
                                            }
                                            if(!(ABHitBuildFlags & (1UL << noABHitInd2))){
                                                addbackE[numABHitsBuilt] += sortedEvt.noABHit[noABHitInd2].energy;
                                            }
                                            //flag hits
                                            ABHitBuildFlags |= (1UL << noABHitInd);
                                            ABHitBuildFlags |= (1UL << noABHitInd2);
                                            ABHitMapping[noABHitInd] = numABHitsBuilt;
                                            ABHitMapping[noABHitInd2] = numABHitsBuilt;
                                            abHitBuilt = 1;
                                        }else{
                                            //first hit not flagged, second hit flagged
                                            //add first hit to second hit's addback data
                                            addbackE[ABHitMapping[noABHitInd2]] += sortedEvt.noABHit[noABHitInd].energy;
                                            //flag hit
                                            ABHitBuildFlags |= (1UL << noABHitInd);
                                            ABHitMapping[noABHitInd] = ABHitMapping[noABHitInd2];
                                        }
                                    }else{
                                        //first hit flagged
                                        if(!(ABHitBuildFlags & (1UL << noABHitInd2))){
                                            //first hit flagged, second not flagged
                                            //add second hit to first hit's addback data
                                            addbackE[ABHitMapping[noABHitInd]] += sortedEvt.noABHit[noABHitInd2].energy;
                                            //flag hit
                                            ABHitBuildFlags |= (1UL << noABHitInd2);
                                            ABHitMapping[noABHitInd2] = ABHitMapping[noABHitInd];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(abHitBuilt != 0){
                    //addback hit was just built
                    numABHitsBuilt++;
                    if(numABHitsBuilt>=64){
                        break;
                    }
                }
            }
        }

        //Check non-addback hits
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            //check coincidences with non-addback hits (addback vs. non-addback was already checked earlier)
            for(uint8_t noABHitInd2=(noABHitInd+1); noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                //first check energy condition
                double sumE = sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy;
                if((sumE >= eLow)&&(sumE<= eHigh)){
                    double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                    //printf("tDiff: %0.3f\n",tDiff);
                    //check if the first hit has already been used
                    uint8_t coincPrevFilled = 0;
                    if(specFillFlags[0] & (1UL << noABHitInd)){
                        coincPrevFilled = 1;
                        break;
                    }
                    if(coincPrevFilled == 0){
                        //check if the second hit has already been used
                        uint8_t coincPrevFilled2 = 0;
                        if(specFillFlags[0] & (1UL << noABHitInd2)){
                            coincPrevFilled2 = 1;
                            break;
                        }
                        if(coincPrevFilled2 == 0){
                            if(tDiff <= COINC_TIMING_GATE){ //timing condition
                                int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                }
                                specFillFlags[0] |= (1UL << noABHitInd); //flag that spectrum was filled with this hit
                                eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                }
                                specFillFlags[0] |= (1UL << noABHitInd2); //flag that spectrum was filled with this hit
                            }
                        }
                    }
                    //fill 180 degree spectra
                    if(!(specFillFlags[1] & (1UL << noABHitInd))){
                        if(tDiff <= COINC_TIMING_GATE){ //timing condition
                            for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                                if(!(specFillFlags[1] & (1UL << noABHitInd3))){
                                    if((noABHitInd3 == noABHitInd)||(noABHitInd3 == noABHitInd2)){
                                        //first hit is a member of either hit already in the sum gate
                                        for(int noABHitInd4 = noABHitInd3+1; noABHitInd4 < sortedEvt.header.numNoABHits; noABHitInd4++){
                                            if(!(specFillFlags[1] & (1UL << noABHitInd4))){
                                                if((noABHitInd4 != noABHitInd)&&(noABHitInd4 != noABHitInd2)){
                                                    //second hit not a member of either hit already in the sum gate
                                                    if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core][sortedEvt.noABHit[noABHitInd4].core] != 0){
                                                        double eGamma = addbackE[ABHitMapping[noABHitInd3]]/keVPerBin;
                                                        double eGamma2 = 0;
                                                        double tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                                        if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                            //opposing hit was part of an addback hit
                                                            eGamma2 = addbackE[ABHitMapping[noABHitInd4]]/keVPerBin;
                                                            if(noABHitInd3 == noABHitInd2){
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - noABHitTime(&sortedEvt,noABHitInd2));
                                                            }else{
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - noABHitTime(&sortedEvt,noABHitInd));
                                                            }
                                                        }else{
                                                            //opposing hit was a single non-addback hit
                                                            eGamma2 = sortedEvt.noABHit[noABHitInd4].energy/keVPerBin;
                                                            if(noABHitInd3 == noABHitInd2){
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - noABHitTime(&sortedEvt,noABHitInd2));
                                                            }else{
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - noABHitTime(&sortedEvt,noABHitInd));
                                                            }
                                                        }
                                                        if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                            //set flags
                                                            specFillFlags[1] |= (1UL << noABHitInd3);
                                                            if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(ABHitMapping[noABHitIndS] == ABHitMapping[noABHitInd4]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }else{
                                                                specFillFlags[1] |= (1UL << noABHitInd4);
                                                            }
                                                            int eGammaSum = (int)(eGamma + eGamma2);
                                                            if(eGammaSum>=0 && eGammaSum<S32K){
                                                                mcaOut[2][eGammaSum]++;
                                                            }
                                                            if(eGamma>=0 && eGamma<S32K){
                                                                mcaOut[1][(int)eGamma]++;
                                                            }
                                                            if(eGamma2>=0 && eGamma2<S32K){
                                                                mcaOut[1][(int)eGamma2]++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
        }

        //for each gate addback hit that falls in the energy gate,
        //build coincident addback hits using the projection addback radius, and
        //add those to the output spectrum
        for(uint8_t ABHitInd=0; ABHitInd < numABHitsBuilt; ABHitInd++){
            //check coincidences with non-addback hits
            for(uint8_t noABHitInd2=0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                //first check energy condition
                double sumE = addbackE[ABHitInd] + sortedEvt.noABHit[noABHitInd2].energy;
                if((sumE >= eLow)&&(sumE<= eHigh)){
                    double tDiff = fabs(addbackT[ABHitInd] - noABHitTime(&sortedEvt,noABHitInd2));
                    //printf("tDiff: %0.3f\n",tDiff);
                    //check if the first hit has already been used
                    uint8_t coincPrevFilled = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                        if(ABHitMapping[noABHitInd] == ABHitInd){
                            if(specFillFlags[0] & (1UL << noABHitInd)){
                                coincPrevFilled = 1;
                                break;
                            }
                        }
                    }
                    if(coincPrevFilled == 0){
                        if(ABHitMapping[noABHitInd2] != ABHitInd){
                            //check if the second hit has already been used
                            uint8_t coincPrevFilled2 = 0;
                            if(specFillFlags[0] & (1UL << noABHitInd2)){
                                coincPrevFilled2 = 1;
                                break;
                            }
                            if(coincPrevFilled2 == 0){
                                if(tDiff <= COINC_TIMING_GATE){ //timing condition
                                    int eGamma = (int)(addbackE[ABHitInd]/keVPerBin);
                                    if(eGamma>=0 && eGamma<S32K){
                                        mcaOut[0][eGamma]++;
                                    }
                                    //flag that spectrum was filled with this hit
                                    for(int noABHitIndB = 0; noABHitIndB < sortedEvt.header.numNoABHits; noABHitIndB++){
                                        if(ABHitMapping[noABHitIndB] == ABHitInd){
                                            specFillFlags[0] |= (1UL << noABHitIndB);
                                        }
                                    }
                                    eGamma = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                    if(eGamma>=0 && eGamma<S32K){
                                        mcaOut[0][eGamma]++;
                                    }
                                    specFillFlags[0] |= (1UL << noABHitInd2); //flag that spectrum was filled with this hit
                                }
                            }
                        }
                    }
                    //fill 180 degree spectra
                    uint8_t coincPrevFilledS = 0;
                    for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                        if(ABHitMapping[noABHitInd3] == ABHitInd){
                            if(specFillFlags[1] & (1UL << noABHitInd3)){
                                coincPrevFilledS = 1;
                                break;
                            }
                        }
                    }
                    if(coincPrevFilledS == 0){
                        if(tDiff <= COINC_TIMING_GATE){ //timing condition
                            for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                                if(!(specFillFlags[1] & (1UL << noABHitInd3))){
                                    if((ABHitMapping[noABHitInd3] == ABHitInd)||(noABHitInd3 == noABHitInd2)){
                                        //first hit is a member of either hit already in the sum gate
                                        for(int noABHitInd4 = noABHitInd3+1; noABHitInd4 < sortedEvt.header.numNoABHits; noABHitInd4++){
                                            if(!(specFillFlags[1] & (1UL << noABHitInd4))){
                                                if((ABHitMapping[noABHitInd4] != ABHitInd)&&(noABHitInd4 != noABHitInd2)){
                                                    //second hit not a member of either hit already in the sum gate
                                                    if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core][sortedEvt.noABHit[noABHitInd4].core] != 0){
                                                        double eGamma = addbackE[ABHitMapping[noABHitInd3]]/keVPerBin;
                                                        double eGamma2 = 0;
                                                        double tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                                        if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                            //opposing hit was part of a different addback hit
                                                            eGamma2 = addbackE[ABHitMapping[noABHitInd4]]/keVPerBin;
                                                            if(noABHitInd3 == noABHitInd2){
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - noABHitTime(&sortedEvt,noABHitInd2));
                                                            }else{
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - addbackT[ABHitInd]);
                                                            }
                                                        }else{
                                                            //opposing hit was a single non-addback hit
                                                            eGamma2 = sortedEvt.noABHit[noABHitInd4].energy/keVPerBin;
                                                            if(noABHitInd3 == noABHitInd2){
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - noABHitTime(&sortedEvt,noABHitInd2));
                                                            }else{
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - addbackT[ABHitInd]);
                                                            }
                                                        }
                                                        if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                            //set flags
                                                            if(noABHitInd3 == noABHitInd2){
                                                                specFillFlags[1] |= (1UL << noABHitInd3);
                                                            }else{
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(ABHitMapping[noABHitIndS] == ABHitMapping[noABHitInd3]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }
                                                            if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(ABHitMapping[noABHitIndS] == ABHitMapping[noABHitInd4]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }else{
                                                                specFillFlags[1] |= (1UL << noABHitInd4);
                                                            }
                                                            int eGammaSum = (int)(eGamma + eGamma2);
                                                            if(eGammaSum>=0 && eGammaSum<S32K){
                                                                mcaOut[2][eGammaSum]++;
                                                            }
                                                            if(eGamma>=0 && eGamma<S32K){
                                                                mcaOut[1][(int)eGamma]++;
                                                            }
                                                            if(eGamma2>=0 && eGamma2<S32K){
                                                                mcaOut[1][(int)eGamma2]++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //check coincidences with other addback hits
            for(uint8_t ABHitInd2=(ABHitInd+1); ABHitInd2 < numABHitsBuilt; ABHitInd2++){
                //first check energy condition
                double sumE = addbackE[ABHitInd] + addbackE[ABHitInd2];
                if((sumE >= eLow)&&(sumE<= eHigh)){
                    double tDiff = fabs(addbackT[ABHitInd] - addbackT[ABHitInd2]);
                    //printf("tDiff: %0.3f\n",tDiff);
                    //check if the first hit has already been used
                    uint8_t coincPrevFilled = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                        if(ABHitMapping[noABHitInd] == ABHitInd){
                            if(specFillFlags[0] & (1UL << noABHitInd)){
                                coincPrevFilled = 1;
                                break;
                            }
                        }
                    }
                    if(coincPrevFilled == 0){
                    
                        //check if the second hit has already been used
                        uint8_t coincPrevFilled2 = 0;
                        for(int noABHitIndB = 0; noABHitIndB < sortedEvt.header.numNoABHits; noABHitIndB++){
                            if(ABHitMapping[noABHitIndB] == ABHitInd2){
                                if(specFillFlags[0] & (1UL << noABHitIndB)){
                                    coincPrevFilled2 = 1;
                                    break;
                                }
                            }
                        }
                        if(coincPrevFilled2 == 0){
                            if(tDiff <= COINC_TIMING_GATE){ //timing condition
                                int eGamma = (int)(addbackE[ABHitInd]/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                }
                                //flag that spectrum was filled with this hit
                                for(int noABHitIndB = 0; noABHitIndB < sortedEvt.header.numNoABHits; noABHitIndB++){
                                    if(ABHitMapping[noABHitIndB] == ABHitInd){
                                        specFillFlags[0] |= (1UL << noABHitIndB);
                                    }
                                }
                                eGamma = (int)(addbackE[ABHitInd2]/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                }
                                //flag that spectrum was filled with this hit
                                for(int noABHitIndB = 0; noABHitIndB < sortedEvt.header.numNoABHits; noABHitIndB++){
                                    if(ABHitMapping[noABHitIndB] == ABHitInd2){
                                        specFillFlags[0] |= (1UL << noABHitIndB);
                                    }
                                }
                                
                            }
                        }
                    }
                    //fill 180 degree spectra
                    uint8_t coincPrevFilledS = 0;
                    for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                        if(ABHitMapping[noABHitInd3] == ABHitInd){
                            if(specFillFlags[1] & (1UL << noABHitInd3)){
                                coincPrevFilledS = 1;
                                break;
                            }
                        }
                    }
                    if(coincPrevFilledS == 0){
                        if(tDiff <= COINC_TIMING_GATE){ //timing condition
                            for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                                if(!(specFillFlags[1] & (1UL << noABHitInd3))){
                                    if((ABHitMapping[noABHitInd3] == ABHitInd)||(ABHitMapping[noABHitInd3] == ABHitInd2)){
                                        //first hit is a member of either hit already in the sum gate
                                        for(int noABHitInd4 = noABHitInd3+1; noABHitInd4 < sortedEvt.header.numNoABHits; noABHitInd4++){
                                            if(!(specFillFlags[1] & (1UL << noABHitInd4))){
                                                if((ABHitMapping[noABHitInd4] != ABHitInd)&&(ABHitMapping[noABHitInd4] != ABHitInd2)){
                                                    //second hit not a member of either hit already in the sum gate
                                                    if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core][sortedEvt.noABHit[noABHitInd4].core] != 0){
                                                        double eGamma = addbackE[ABHitMapping[noABHitInd3]]/keVPerBin;
                                                        double eGamma2 = 0;
                                                        double tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                                        if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                            //opposing hit was part of a different addback hit
                                                            eGamma2 = addbackE[ABHitMapping[noABHitInd4]]/keVPerBin;
                                                            if(ABHitMapping[noABHitInd3] == ABHitInd){
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - addbackT[ABHitInd]);
                                                            }else{
                                                                tDiff = fabs(addbackT[ABHitMapping[noABHitInd4]] - addbackT[ABHitInd2]);
                                                            }
                                                        }else{
                                                            //opposing hit was a single non-addback hit
                                                            eGamma2 = sortedEvt.noABHit[noABHitInd4].energy/keVPerBin;
                                                            if(ABHitMapping[noABHitInd3] == ABHitInd){
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - addbackT[ABHitInd]);
                                                            }else{
                                                                tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd4) - addbackT[ABHitInd2]);
                                                            }
                                                        }
                                                        if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                            //set flags
                                                            for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                if(ABHitMapping[noABHitIndS] == ABHitMapping[noABHitInd3]){
                                                                    specFillFlags[1] |= (1UL << noABHitIndS);
                                                                }
                                                            }
                                                            if(ABHitBuildFlags & (1UL << noABHitInd4)){
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(ABHitMapping[noABHitIndS] == ABHitMapping[noABHitInd4]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }else{
                                                                specFillFlags[1] |= (1UL << noABHitInd4);
                                                            }
                                                            int eGammaSum = (int)(eGamma + eGamma2);
                                                            if(eGammaSum>=0 && eGammaSum<S32K){
                                                                mcaOut[2][eGammaSum]++;
                                                            }
                                                            if(eGamma>=0 && eGamma<S32K){
                                                                mcaOut[1][(int)eGamma]++;
                                                            }
                                                            if(eGamma2>=0 && eGamma2<S32K){
                                                                mcaOut[1][(int)eGamma2]++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
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

  EEGamma_modAB_sumEgate_mca_SMOL *mysort = new EEGamma_modAB_sumEgate_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh, ABRad;
  printf("Starting EEGamma_modAB_sumEgate_mca_SMOL\n");

  if((argc != 6)&&(argc != 7)){
    cout << "Generates dmca gamma-gated spectra." << endl;
    cout << "Arguments: EEGamma_modAB_sumEgate_mca_SMOL smol_file EGateLow EGateHigh ABRadius output_dmca_file keV_per_bin" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *ABRadius* is in mm. A value of zero corresponds to no addback." << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    eLow = atof(argv[2]);
    eHigh = atof(argv[3]);
    ABRad = atof(argv[4]);
    outfile = argv[5];
    if(argc > 6){
      keVPerBin = atof(argv[6]);
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
  if(ABRad < 0.0){
    cout << "ERROR: addback radius must be 0 or larger!" << endl;
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
    printf("SMOL tree: %s\nEnergy gate: [%0.2f %0.2f]\nGate addback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, ABRad, outfile, keVPerBin);
    numSepEvts += mysort->SortData(sfile, eLow, eHigh, ABRad, keVPerBin);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nEnergy gate: [%0.2f %0.2f]\nGate addback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, ABRad, outfile, keVPerBin);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, eLow, eHigh, ABRad, keVPerBin);
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
