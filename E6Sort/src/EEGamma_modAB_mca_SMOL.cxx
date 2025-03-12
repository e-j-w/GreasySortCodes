//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_modAB_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_modAB_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

uint8_t gateABHitMapping[64], projABHitMapping[64]; //arrays specifying which hits correspond to which addback hits
double gateAddbackE[64], projAddbackE[64];
double gateAddbackT[64], projAddbackT[64];

void EEGamma_modAB_mca_SMOL::WriteData(const char* outName){

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

uint64_t EEGamma_modAB_mca_SMOL::SortData(const char *sfile, const double eLow, const double eHigh, const double gateABRad, const double projABRad, const double keVPerBin){

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
                        if(getGeHitDistance(i,0,k,0,1) < projABRad){
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
        uint64_t gateABHitBuildFlags = 0;
        uint8_t numGateABHitsBuilt = 0;
        memset(gateABHitMapping,0,sizeof(gateABHitMapping));
        memset(gateAddbackE,0,sizeof(gateAddbackE));
        memset(gateAddbackT,0,sizeof(gateAddbackT));

        //build initial addback hits, using the gate addback radius,
        //and using gateABHitMapping to track which hits are grouped together
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
        if(noABHitInd < 64){
            uint8_t abHitBuilt = 0;
            for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                if(noABHitInd2 < 64){
                    if(noABHitInd2 != noABHitInd){
                        //check if hits are in neighbouring crystals
                        if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < gateABRad){ //FORWARD POSITION (11 cm)
                            double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                            if(tDiff <= ADDBACK_TIMING_GATE){ //timing condition
                                if(!(gateABHitBuildFlags & (1UL << noABHitInd))){
                                    //first hit not yet flagged
                                    if(!(gateABHitBuildFlags & (1UL << noABHitInd2))){
                                        //neither hit flagged
                                        //build initial addback hit
                                        gateAddbackT[numGateABHitsBuilt] = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                                        if(!(gateABHitBuildFlags & (1UL << noABHitInd))){
                                            gateAddbackE[numGateABHitsBuilt] += sortedEvt.noABHit[noABHitInd].energy;
                                        }
                                        if(!(gateABHitBuildFlags & (1UL << noABHitInd2))){
                                            gateAddbackE[numGateABHitsBuilt] += sortedEvt.noABHit[noABHitInd2].energy;
                                        }
                                        //flag hits
                                        gateABHitBuildFlags |= (1UL << noABHitInd);
                                        gateABHitBuildFlags |= (1UL << noABHitInd2);
                                        gateABHitMapping[noABHitInd] = numGateABHitsBuilt;
                                        gateABHitMapping[noABHitInd2] = numGateABHitsBuilt;
                                        abHitBuilt = 1;
                                    }else{
                                        //first hit not flagged, second hit flagged
                                        //add first hit to second hit's addback data
                                        gateAddbackE[gateABHitMapping[noABHitInd2]] += sortedEvt.noABHit[noABHitInd].energy;
                                        //flag hit
                                        gateABHitBuildFlags |= (1UL << noABHitInd);
                                        gateABHitMapping[noABHitInd] = gateABHitMapping[noABHitInd2];
                                    }
                                }else{
                                    //first hit flagged
                                    if(!(gateABHitBuildFlags & (1UL << noABHitInd2))){
                                        //first hit flagged, second not flagged
                                        //add second hit to first hit's addback data
                                        gateAddbackE[gateABHitMapping[noABHitInd]] += sortedEvt.noABHit[noABHitInd2].energy;
                                        //flag hit
                                        gateABHitBuildFlags |= (1UL << noABHitInd2);
                                        gateABHitMapping[noABHitInd2] = gateABHitMapping[noABHitInd];
                                    }
                                }
                            }
                        }
                    }

                }
            }
            if(abHitBuilt != 0){
                //addback hit was just built
                numGateABHitsBuilt++;
                if(numGateABHitsBuilt>=64){
                    break;
                }
            }
        }
        }
        //printf("num original hits: %u\n",sortedEvt.header.numNoABHits);
        //printf("  numGateHits: %u\n",numGateABHitsBuilt);

        //for each gate addback hit that falls in the energy gate,
        //build coincident addback hits using the projection addback radius, and
        //add those to the output spectrum
        uint8_t specFilled = 0;
        for(uint8_t gateABHitInd=0; gateABHitInd < numGateABHitsBuilt; gateABHitInd++){
            //evaluate whether or not any of the individual hits used to construct the gate
            //addback hit have been flagged as already used in the final spectrum
            if((gateAddbackE[gateABHitInd] >= eLow)&&(gateAddbackE[gateABHitInd] <= eHigh)){
                //satisfactory gate hit found, flag its original hits as in use
                uint64_t usedGateHitBuildFlags = 0;
                for(uint8_t i=0;i<sortedEvt.header.numNoABHits;i++){
                    if(i<64){
                        if(gateABHitBuildFlags & (1UL << i)){
                            if(gateABHitMapping[i]==gateABHitInd){
                                usedGateHitBuildFlags |= (1UL << i);
                            }
                        }
                    }
                }
                //now, look at all other hits, and constuct projection addback hits using them
                uint64_t projABHitBuildFlags = 0;
                uint8_t numProjABHitsBuilt = 0;
                memset(projABHitMapping,0,sizeof(projABHitMapping));
                memset(projAddbackE,0,sizeof(projAddbackE));
                memset(projAddbackT,0,sizeof(projAddbackT));
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    if(noABHitInd < 64){
                        if(!(usedGateHitBuildFlags & (1UL << noABHitInd))){
                            uint8_t abHitBuilt = 0;
                            for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                if(noABHitInd2 < 64){
                                    if((!(usedGateHitBuildFlags & (1UL << noABHitInd2)))&&(noABHitInd2 != noABHitInd)){
                                        //check if hits are in neighbouring crystals
                                        if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < projABRad){ //FORWARD POSITION (11 cm)
                                            double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                            if(tDiff <= ADDBACK_TIMING_GATE){ //timing condition
                                                if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                    //first hit not yet flagged
                                                    if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                        //neither hit flagged
                                                        //build initial addback hit
                                                        projAddbackT[numProjABHitsBuilt] = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                            projAddbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd].energy;
                                                        }
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                            projAddbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd2].energy;
                                                        }
                                                        //flag hits
                                                        projABHitBuildFlags |= (1UL << noABHitInd);
                                                        projABHitBuildFlags |= (1UL << noABHitInd2);
                                                        projABHitMapping[noABHitInd] = numProjABHitsBuilt;
                                                        projABHitMapping[noABHitInd2] = numProjABHitsBuilt;
                                                        abHitBuilt = 1;
                                                    }else{
                                                        //first hit not flagged, second hit flagged
                                                        //add first hit to second hit's addback data
                                                        projAddbackE[projABHitMapping[noABHitInd2]] += sortedEvt.noABHit[noABHitInd].energy;
                                                        //flag hit
                                                        projABHitBuildFlags |= (1UL << noABHitInd);
                                                        projABHitMapping[noABHitInd] = projABHitMapping[noABHitInd2];
                                                    }
                                                }else{
                                                    //first hit flagged
                                                    if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                        //first hit flagged, second not flagged
                                                        //add second hit to first hit's addback data
                                                        projAddbackE[projABHitMapping[noABHitInd]] += sortedEvt.noABHit[noABHitInd2].energy;
                                                        //flag hit
                                                        projABHitBuildFlags |= (1UL << noABHitInd2);
                                                        projABHitMapping[noABHitInd2] = projABHitMapping[noABHitInd];
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if(abHitBuilt != 0){
                                //addback hit was just built
                                numProjABHitsBuilt++;
                                if(numProjABHitsBuilt>=64){
                                    break;
                                }
                            }
                        }
                    }
                }
                //projection addback hits built, now add them to the spectrum
                //printf("  numProjHits: %u\n",numProjABHitsBuilt);
                for(uint8_t projABHitInd=0; projABHitInd < numProjABHitsBuilt; projABHitInd++){
                    double tDiff = fabs(gateAddbackT[gateABHitInd] - projAddbackT[projABHitInd]);
                    if(tDiff<= COINC_TIMING_GATE){ //timing condition
                        int eGamma = (int)(projAddbackE[projABHitInd]/keVPerBin);
                        if(eGamma>=0 && eGamma<S32K){
                            mcaOut[0][eGamma]++;
                        }
                        //check angles of original projection hits, and fill 180 degree spectra
                        uint8_t oppHitFound = 0;
                        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                            if(projABHitMapping[noABHitInd] == projABHitInd){
                                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                    if(((!(projABHitBuildFlags & (1UL << noABHitInd2))) || (projABHitMapping[noABHitInd2] != projABHitInd))&&((!(usedGateHitBuildFlags & (1UL << noABHitInd2))) || (gateABHitMapping[noABHitInd2] != gateABHitInd))){
                                        if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                                            int eGamma2 = 0;
                                            tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                            if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                //opposing hit was part of a different projection hit
                                                eGamma2 = (int)(projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin);
                                                tDiff = fabs(projAddbackT[projABHitMapping[noABHitInd2]] - projAddbackT[projABHitInd]);
                                            }else{
                                                //opposing hit was a single non-addback hit
                                                eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                                tDiff = fabs(sortedEvt.noABHit[noABHitInd2].timeOffsetNs - projAddbackT[projABHitInd]);
                                            }
                                            if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                int eGammaSum = eGamma + eGamma2;
                                                if(eGammaSum>=0 && eGammaSum<S32K){
                                                    mcaOut[2][eGammaSum]++;
                                                }
                                                if(eGamma>=0 && eGamma<S32K){
                                                    mcaOut[1][eGamma]++;
                                                }
                                                if(eGamma2>=0 && eGamma2<S32K){
                                                    mcaOut[1][eGamma2]++;
                                                }
                                                oppHitFound = 1;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                            if(oppHitFound){
                                break;
                            }
                        }
                    }else{
                        //time random
                        int eGamma = (int)(projAddbackE[projABHitInd]/keVPerBin);
                        if(eGamma>=0 && eGamma<S32K){
                            mcaOut[3][eGamma]++;
                        }
                    }
                }
                //add in any remaining hits that weren't addback'd
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    if(!(usedGateHitBuildFlags & (1UL << noABHitInd))){
                        if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                            double tDiff = fabs(gateAddbackT[gateABHitInd] - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                            if(tDiff<= COINC_TIMING_GATE){ //timing condition
                                int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                }
                                //check angles of original projection hits, and fill 180 degree spectra
                                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                    if((!(usedGateHitBuildFlags & (1UL << noABHitInd2))) || (gateABHitMapping[noABHitInd2] != gateABHitInd)){
                                        if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                                            int eGamma2 = 0;
                                            tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                            if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                //opposing hit was part of a different projection hit
                                                eGamma2 = (int)(projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin);
                                                tDiff = fabs(projAddbackT[projABHitMapping[noABHitInd2]] - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                                            }else{
                                                //opposing hit was a single non-addback hit
                                                eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                                tDiff = fabs(sortedEvt.noABHit[noABHitInd2].timeOffsetNs - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                                            }
                                            if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                int eGammaSum = eGamma + eGamma2;
                                                if(eGammaSum>=0 && eGammaSum<S32K){
                                                    mcaOut[2][eGammaSum]++;
                                                }
                                                if(eGamma>=0 && eGamma<S32K){
                                                    mcaOut[1][eGamma]++;
                                                }
                                                if(eGamma2>=0 && eGamma2<S32K){
                                                    mcaOut[1][eGamma2]++;
                                                }
                                                break;
                                            }
                                        }
                                    }
                                }
                            }else{
                                //time random
                                int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[3][eGamma]++;
                                }
                            }
                        }
                    }
                }
                specFilled = 1;
                break; //don't need to consider any other gammas in the energy gate, already added all coincident gammas to the spectrum
            }
        }

        if(specFilled == 0){
            //no addback hits of appropriate energy (maybe none at all), check non-addback hits
            for(int gateNoABHitInd = 0; gateNoABHitInd < sortedEvt.header.numNoABHits; gateNoABHitInd++){
                if((sortedEvt.noABHit[gateNoABHitInd].energy >= eLow)&&(sortedEvt.noABHit[gateNoABHitInd].energy <= eHigh)){
                    //now, look at all other hits, and constuct projection addback hits using them
                    uint64_t projABHitBuildFlags = 0;
                    uint8_t numProjABHitsBuilt = 0;
                    memset(projABHitMapping,0,sizeof(projABHitMapping));
                    memset(projAddbackE,0,sizeof(projAddbackE));
                    memset(projAddbackT,0,sizeof(projAddbackT));
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                        if(noABHitInd < 64){
                            if(noABHitInd != gateNoABHitInd){
                                uint8_t abHitBuilt = 0;
                                for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                    if(noABHitInd2 < 64){
                                        if((noABHitInd2 != gateNoABHitInd)&&(noABHitInd2 != noABHitInd)){
                                            //check if hits are in neighbouring crystals
                                            if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < projABRad){ //FORWARD POSITION (11 cm)
                                                double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                                if(tDiff <= ADDBACK_TIMING_GATE){ //timing condition
                                                    if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                        //first hit not yet flagged
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                            //neither hit flagged
                                                            //build initial addback hit
                                                            projAddbackT[numProjABHitsBuilt] = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                                                            if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                                projAddbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd].energy;
                                                            }
                                                            if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                                projAddbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd2].energy;
                                                            }
                                                            //flag hits
                                                            projABHitBuildFlags |= (1UL << noABHitInd);
                                                            projABHitBuildFlags |= (1UL << noABHitInd2);
                                                            projABHitMapping[noABHitInd] = numProjABHitsBuilt;
                                                            projABHitMapping[noABHitInd2] = numProjABHitsBuilt;
                                                            abHitBuilt = 1;
                                                        }else{
                                                            //first hit not flagged, second hit flagged
                                                            //add first hit to second hit's addback data
                                                            projAddbackE[projABHitMapping[noABHitInd2]] += sortedEvt.noABHit[noABHitInd].energy;
                                                            //flag hit
                                                            projABHitBuildFlags |= (1UL << noABHitInd);
                                                            projABHitMapping[noABHitInd] = projABHitMapping[noABHitInd2];
                                                        }
                                                    }else{
                                                        //first hit flagged
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                            //first hit flagged, second not flagged
                                                            //add second hit to first hit's addback data
                                                            projAddbackE[projABHitMapping[noABHitInd]] += sortedEvt.noABHit[noABHitInd2].energy;
                                                            //flag hit
                                                            projABHitBuildFlags |= (1UL << noABHitInd2);
                                                            projABHitMapping[noABHitInd2] = projABHitMapping[noABHitInd];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                    }
                                }
                                if(abHitBuilt != 0){
                                    //addback hit was just built
                                    numProjABHitsBuilt++;
                                    if(numProjABHitsBuilt>=64){
                                        break;
                                    }
                                }
                            }
                            
                        }
                        
                    }
                    //projection addback hits built, now add them to the spectrum
                    //printf("  numProjHits: %u\n",numProjABHitsBuilt);
                    for(uint8_t projABHitInd=0; projABHitInd < numProjABHitsBuilt; projABHitInd++){
                        double tDiff = fabs(sortedEvt.noABHit[gateNoABHitInd].timeOffsetNs - projAddbackT[projABHitInd]);
                        if(tDiff<= COINC_TIMING_GATE){ //timing condition
                            int eGamma = (int)(projAddbackE[projABHitInd]/keVPerBin);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[0][eGamma]++;
                            }
                            //check angles of original projection hits, and fill 180 degree spectra
                            uint8_t oppHitFound = 0;
                            for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                                if(projABHitMapping[noABHitInd] == projABHitInd){
                                    for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                        if(((!(projABHitBuildFlags & (1UL << noABHitInd2))) || (projABHitMapping[noABHitInd2] != projABHitInd))&&(noABHitInd2 != gateNoABHitInd)){
                                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                                                int eGamma2 = 0;
                                                tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                                if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                    //opposing hit was part of a different projection hit
                                                    eGamma2 = (int)(projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin);
                                                    tDiff = fabs(projAddbackT[projABHitMapping[noABHitInd2]] - projAddbackT[projABHitInd]);
                                                }else{
                                                    //opposing hit was a single non-addback hit
                                                    eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                                    tDiff = fabs(sortedEvt.noABHit[noABHitInd2].timeOffsetNs - projAddbackT[projABHitInd]);
                                                }
                                                if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                    int eGammaSum = eGamma + eGamma2;
                                                    if(eGammaSum>=0 && eGammaSum<S32K){
                                                        mcaOut[2][eGammaSum]++;
                                                    }
                                                    if(eGamma>=0 && eGamma<S32K){
                                                        mcaOut[1][eGamma]++;
                                                    }
                                                    if(eGamma2>=0 && eGamma2<S32K){
                                                        mcaOut[1][eGamma2]++;
                                                    }
                                                    oppHitFound = 1;
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                if(oppHitFound){
                                    break;
                                }
                            }
                        }else{
                            //time random
                            int eGamma = (int)(projAddbackE[projABHitInd]/keVPerBin);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[3][eGamma]++;
                            }
                        }
                    }
                    //add in any remaining hits that weren't addback'd
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                        if(noABHitInd != gateNoABHitInd){
                            if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                double tDiff = fabs(sortedEvt.noABHit[gateNoABHitInd].timeOffsetNs - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                                if(tDiff<= COINC_TIMING_GATE){ //timing condition
                                    int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                    if(eGamma>=0 && eGamma<S32K){
                                        mcaOut[0][eGamma]++;
                                    }
                                    //check angles of original projection hits, and fill 180 degree spectra
                                    for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                        if((noABHitInd2 != noABHitInd)&&(noABHitInd2 != gateNoABHitInd)){
                                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                                                int eGamma2 = 0;
                                                tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                                if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                    //opposing hit was part of a different projection hit
                                                    eGamma2 = (int)(projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin);
                                                    tDiff = fabs(projAddbackT[projABHitMapping[noABHitInd2]] - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                                                }else{
                                                    //opposing hit was a single non-addback hit
                                                    eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                                    tDiff = fabs(sortedEvt.noABHit[noABHitInd2].timeOffsetNs - sortedEvt.noABHit[noABHitInd].timeOffsetNs);
                                                }
                                                if(tDiff<= SUM_TIMING_GATE){ //timing condition
                                                    int eGammaSum = eGamma + eGamma2;
                                                    if(eGammaSum>=0 && eGammaSum<S32K){
                                                        mcaOut[2][eGammaSum]++;
                                                    }
                                                    if(eGamma>=0 && eGamma<S32K){
                                                        mcaOut[1][eGamma]++;
                                                    }
                                                    if(eGamma2>=0 && eGamma2<S32K){
                                                        mcaOut[1][eGamma2]++;
                                                    }
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }else{
                                    //time random
                                    int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                    if(eGamma>=0 && eGamma<S32K){
                                        mcaOut[3][eGamma]++;
                                    }
                                }
                            }
                        }
                    }
                    specFilled = 1;
                    break; //don't need to consider any other gammas in the energy gate, already added all coincident gammas to the spectrum
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

  EEGamma_modAB_mca_SMOL *mysort = new EEGamma_modAB_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh, gateABRad, projABRad;
  printf("Starting EEGamma_modAB_mca_SMOL\n");

  if((argc != 7)&&(argc != 8)){
    cout << "Generates dmca gamma-gated spectra." << endl;
    cout << "Arguments: EEGamma_modAB_mca_SMOL smol_file EGateLow EGateHigh gateABRadius projABRadius output_dmca_file keV_per_bin" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *gateABRadius* and *projABRadius* are in mm. Values of zero correspond to no addback" << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    eLow = atof(argv[2]);
    eHigh = atof(argv[3]);
    gateABRad = atof(argv[4]);
    projABRad = atof(argv[5]);
    outfile = argv[6];
    if(argc > 7){
      keVPerBin = atof(argv[7]);
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
  if((gateABRad < 0.0)||(projABRad < 0.0)){
    cout << "ERROR: addback radii must be 0 or larger!" << endl;
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
    printf("SMOL tree: %s\nEnergy gate: [%0.2f %0.2f]\nGate addback radius: %0.2f mm\nProjection addback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, gateABRad, projABRad, outfile, keVPerBin);
    numSepEvts += mysort->SortData(sfile, eLow, eHigh, gateABRad, projABRad, keVPerBin);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nEnergy gate: [%0.2f %0.2f]\nGate addback radius: %0.2f mm\nProjection addback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, gateABRad, projABRad, outfile, keVPerBin);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, eLow, eHigh, gateABRad, projABRad, keVPerBin);
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
