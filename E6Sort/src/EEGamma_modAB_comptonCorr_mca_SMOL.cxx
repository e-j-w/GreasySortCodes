//Generates addback gamma ray spectra, constructing hits on-the-fly
//with user-defined addback windows for gate and projection hits.
//timing windows are defined in common.h

//Double filling of spectra is avoided using the specFillFlags bit-pattern.
//It is possible to have the gate hit be one of the already filled hits,
//which is necessary if multiple hits fall into the energy gate.

#define EEGamma_modAB_comptonCorr_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_modAB_comptonCorr_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

uint8_t gateABHitMapping[64], projABHitMapping[64]; //arrays specifying which hits correspond to which addback hits
double gateAddbackE[64], projAddbackE[64];
double gateAddbackT[64], projAddbackT[64];

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

uint64_t SortData(const char *sfile, const double eLow, const double eHigh, const double eShift, const double gateABRad, const double projABRad, const uint8_t sameCloverOnly, const uint8_t forwardPos, const double keVPerBin, const uint8_t discardPileup){

    FILE *inp = fopen(sfile, "rb");
    printf("\nFile %s opened\n", sfile);
    
    uint64_t sentries = 0U;
    uint64_t pileupCtrs[16];
    fread(&sentries,sizeof(uint64_t),1,inp);
    uint64_t smolVersion = (uint64_t)(sentries >> 48);
    if(smolVersion > 0){
        fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
        //printf("\nNumber of hits of each pileup type:\n");
        uint64_t totalHits = 0;
        for(uint8_t i=0; i<16; i++){
        //printf("Pileup type %2u: %Lu\n",i,pileupCtrs[i]);
        totalHits += pileupCtrs[i];
        }
        //printf("Total hits:     %Lu\n",totalHits);
        long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
        printf("Fraction of hits with pileup type 1 (no pileup): %Lf\n",frac);
    }
    sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
    sorted_evt sortedEvt;

    //construct 180 degree summing hit map
    memset(hitMap180deg,0,sizeof(hitMap180deg));
    for(uint8_t i=0;i<64;i++){ //first core
        for(uint8_t j=0;j<64;j++){ //coinc core
            if(i!=j){
                if(getGRIFFINVector(i,1).Angle(getGRIFFINVector(j,1))*180.0/PI > 175.0){
                    hitMap180deg[i][j] = 1;
                    continue; //check the next coinc core
                }
                for(uint8_t k=0;k<64;k++){ //other core next to the first core, which could be addback'd with ti
                    if((k!=i)&&(k!=j)){
                        if((!sameCloverOnly) || ((i/4) == (k/4))){
                            if(getGRIFFINHitDistance(i,k,forwardPos) < projABRad){
                                if(getGRIFFINVector(k,1).Angle(getGRIFFINVector(j,1))*180.0/PI > 175.0){
                                    hitMap180deg[i][j] = 1;
                                    break; //check the next coinc core
                                }
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

                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                    continue; //skip pileup hit
                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                    continue; //skip non-pileup hit
                }

                uint8_t abHitBuilt = 0;
                for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                    if(noABHitInd2 < 64){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(noABHitInd2 != noABHitInd){
                            //check if hits are in neighbouring crystals
                            if((!sameCloverOnly) || (((sortedEvt.noABHit[noABHitInd].core & 63U)/4) == ((sortedEvt.noABHit[noABHitInd2].core & 63U)/4))){
                                if(getGRIFFINHitDistance(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos) < gateABRad){ //FORWARD POSITION (11 cm)
                                    double tDiff = (noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                    if(fabs(tDiff) <= ADDBACK_TIMING_GATE){ //timing condition
                                        if(!(gateABHitBuildFlags & (1UL << noABHitInd))){
                                            //first hit not yet flagged
                                            if(!(gateABHitBuildFlags & (1UL << noABHitInd2))){
                                                //neither hit flagged
                                                //build initial addback hit
                                                gateAddbackT[numGateABHitsBuilt] = noABHitTime(&sortedEvt,noABHitInd);
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

        uint64_t specFillFlags[2]; //flags specifying which hits have been used to fill spectra
        memset(specFillFlags,0,sizeof(specFillFlags));

        //check coincidences with non-addback hits first
        //(could look for coincidences with addback hits first, but empirically this approach seems to give lower background)
        for(int gateNoABHitInd = 0; gateNoABHitInd < sortedEvt.header.numNoABHits; gateNoABHitInd++){
            if((sortedEvt.noABHit[gateNoABHitInd].energy >= eLow)&&(sortedEvt.noABHit[gateNoABHitInd].energy <= eHigh)){
                //look at all other hits, and constuct projection addback hits using them
                uint64_t projABHitBuildFlags = 0;
                uint8_t numProjABHitsBuilt = 0;
                memset(projABHitMapping,0,sizeof(projABHitMapping));
                memset(projAddbackE,0,sizeof(projAddbackE));
                memset(projAddbackT,0,sizeof(projAddbackT));
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    if(noABHitInd < 64){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(noABHitInd != gateNoABHitInd){
                            uint8_t abHitBuilt = 0;
                            for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                if(noABHitInd2 < 64){

                                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                        continue; //skip pileup hit
                                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                        continue; //skip non-pileup hit
                                    }

                                    if((noABHitInd2 != gateNoABHitInd)&&(noABHitInd2 != noABHitInd)){
                                        //check if hits are in neighbouring crystals
                                        if((!sameCloverOnly) || (((sortedEvt.noABHit[noABHitInd].core & 63U)/4) == ((sortedEvt.noABHit[noABHitInd2].core & 63U)/4))){
                                            if(getGRIFFINHitDistance(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos) < projABRad){ //FORWARD POSITION (11 cm)
                                                double tDiff = (noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                                if(fabs(tDiff) <= ADDBACK_TIMING_GATE){ //timing condition
                                                    if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                        //first hit not yet flagged
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                            //neither hit flagged
                                                            //build initial addback hit
                                                            projAddbackT[numProjABHitsBuilt] = noABHitTime(&sortedEvt,noABHitInd);
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
                    double tDiff = (noABHitTime(&sortedEvt,gateNoABHitInd) - projAddbackT[projABHitInd]);
                    //printf("tDiff: %0.3f\n",tDiff);
                    //check flags
                    uint8_t hitPrevFilled = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(projABHitMapping[noABHitInd] == projABHitInd){
                            if(specFillFlags[0] & (1UL << noABHitInd)){
                                hitPrevFilled = 1;
                                break;
                            }
                        }
                    }
                    if(hitPrevFilled == 0){
                        if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                            //fill spectrum with hit
                            int eGamma = (int)((projAddbackE[projABHitInd]/keVPerBin) + eShift);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[0][eGamma]++;
                                //mcaOut[6][eGamma]++;
                            }
                            //flag that spectrum was filled with this hit
                            for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                                    continue; //skip pileup hit
                                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                                    continue; //skip non-pileup hit
                                }

                                if(projABHitMapping[noABHitInd] == projABHitInd){
                                    specFillFlags[0] |= (1UL << noABHitInd);
                                }
                            }
                        }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                            //time random
                            int eGamma = (int)((projAddbackE[projABHitInd]/keVPerBin) + eShift);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[3][eGamma]++;
                            }
                        }
                        
                    }
                    //check angles of original projection hits, and fill 180 degree spectra
                    //first, check flags
                    uint8_t hitPrevFilledS = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(projABHitMapping[noABHitInd] == projABHitInd){
                            if(specFillFlags[1] & (1UL << noABHitInd)){
                                hitPrevFilledS = 1;
                                break;
                            }
                        }
                    }
                    if(hitPrevFilledS == 0){
                        if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                            for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                                    continue; //skip pileup hit
                                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                                    continue; //skip non-pileup hit
                                }

                                if(!(specFillFlags[1] & (1UL << noABHitInd))){
                                    if(projABHitMapping[noABHitInd] == projABHitInd){
                                        for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

                                            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                                continue; //skip pileup hit
                                            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                                continue; //skip non-pileup hit
                                            }

                                            if(!(specFillFlags[1] & (1UL << noABHitInd2))){
                                                if(((!(projABHitBuildFlags & (1UL << noABHitInd2))) || (projABHitMapping[noABHitInd2] != projABHitInd))&&(noABHitInd2 != gateNoABHitInd)){
                                                    if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                                        double eGamma = (projAddbackE[projABHitInd]/keVPerBin) + eShift;
                                                        double eGamma2 = 0;
                                                        tDiff = SUM_TIMING_GATE_MAX + 1000.0; //default value, outside the gate
                                                        if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                            //opposing hit was part of a different projection hit
                                                            eGamma2 = (projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin) + eShift;
                                                            tDiff = (projAddbackT[projABHitMapping[noABHitInd2]] - projAddbackT[projABHitInd]);
                                                        }else{
                                                            //opposing hit was a single non-addback hit
                                                            eGamma2 = (sortedEvt.noABHit[noABHitInd2].energy/keVPerBin) + eShift;
                                                            tDiff = (noABHitTime(&sortedEvt,noABHitInd2) - projAddbackT[projABHitInd]);
                                                        }
                                                        if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                                            //set flags
                                                            for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                if(projABHitMapping[noABHitIndS] == projABHitInd){
                                                                    specFillFlags[1] |= (1UL << noABHitIndS);
                                                                }
                                                            }
                                                            if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(projABHitMapping[noABHitIndS] == projABHitMapping[noABHitInd2]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }else{
                                                                specFillFlags[1] |= (1UL << noABHitInd2);
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
                
                //add in any remaining hits that weren't addback'd
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                        continue; //skip pileup hit
                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                        continue; //skip non-pileup hit
                    }

                    if(noABHitInd != gateNoABHitInd){
                        double tDiff = (noABHitTime(&sortedEvt,gateNoABHitInd) - noABHitTime(&sortedEvt,noABHitInd));
                        //printf("tDiff: %0.3f\n",tDiff);
                        if(!(specFillFlags[0] & (1UL << noABHitInd))){
                            if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                                //fill spectrum with hit
                                int eGamma = (int)((sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                    //mcaOut[7][eGamma]++;
                                }
                                //flag that spectrum was filled with this hit
                                specFillFlags[0] |= (1UL << noABHitInd);
                            }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                                //time random
                                int eGamma = (int)((sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[3][eGamma]++;
                                }
                            }
                        }
                        //check angles of original projection hits, and fill 180 degree spectra
                        if(!(specFillFlags[1] & (1UL << noABHitInd))){
                            if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

                                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                        continue; //skip pileup hit
                                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                        continue; //skip non-pileup hit
                                    }

                                    if(!(specFillFlags[1] & (1UL << noABHitInd2))){
                                        if((noABHitInd2 != noABHitInd)&&(noABHitInd2 != gateNoABHitInd)){
                                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                                double eGamma = (sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift;
                                                double eGamma2 = 0;
                                                double tDiff = SUM_TIMING_GATE_MAX + 1000.0; //default value, outside the gate
                                                if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                    //opposing hit was part of a different projection hit
                                                    eGamma2 = (projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin) + eShift;
                                                    tDiff = (projAddbackT[projABHitMapping[noABHitInd2]] - noABHitTime(&sortedEvt,noABHitInd));
                                                }else{
                                                    //opposing hit was a single non-addback hit
                                                    eGamma2 = (sortedEvt.noABHit[noABHitInd2].energy/keVPerBin) + eShift;
                                                    tDiff = (noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
                                                }
                                                if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                                    //set flags
                                                    specFillFlags[1] |= (1UL << noABHitInd);
                                                    if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                        for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                            if(projABHitMapping[noABHitIndS] == projABHitMapping[noABHitInd2]){
                                                                specFillFlags[1] |= (1UL << noABHitIndS);
                                                            }
                                                        }
                                                    }else{
                                                        specFillFlags[1] |= (1UL << noABHitInd2);
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
                //need to continue to later iterations of the loop, in case there are multiple coincident pairs
            }
        }

        //check coincidences with addback hits
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
                //now, look at all other hits, and constuct projection addback hits from them
                //using the projection addback radius
                uint64_t projABHitBuildFlags = 0;
                uint8_t numProjABHitsBuilt = 0;
                memset(projABHitMapping,0,sizeof(projABHitMapping));
                memset(projAddbackE,0,sizeof(projAddbackE));
                memset(projAddbackT,0,sizeof(projAddbackT));
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
                    if(noABHitInd < 64){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(!(usedGateHitBuildFlags & (1UL << noABHitInd))){
                            uint8_t abHitBuilt = 0;
                            for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                                if(noABHitInd2 < 64){

                                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                        continue; //skip pileup hit
                                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                        continue; //skip non-pileup hit
                                    }

                                    if((!(usedGateHitBuildFlags & (1UL << noABHitInd2)))&&(noABHitInd2 != noABHitInd)){
                                        //check if hits are in neighbouring crystals
                                        if((!sameCloverOnly) || (((sortedEvt.noABHit[noABHitInd].core & 63U)/4) == ((sortedEvt.noABHit[noABHitInd2].core & 63U)/4))){
                                            if(getGRIFFINHitDistance(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos) < projABRad){ //FORWARD POSITION (11 cm)
                                                double tDiff = (noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                                if(fabs(tDiff) <= ADDBACK_TIMING_GATE){ //timing condition
                                                    if(!(projABHitBuildFlags & (1UL << noABHitInd))){
                                                        //first hit not yet flagged
                                                        if(!(projABHitBuildFlags & (1UL << noABHitInd2))){
                                                            //neither hit flagged
                                                            //build initial addback hit
                                                            projAddbackT[numProjABHitsBuilt] = noABHitTime(&sortedEvt,noABHitInd);
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
                    double tDiff = (gateAddbackT[gateABHitInd] - projAddbackT[projABHitInd]);
                    //printf("tDiff: %0.3f\n",tDiff);
                    //check flags
                    uint8_t hitPrevFilled = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(projABHitMapping[noABHitInd] == projABHitInd){
                            if(specFillFlags[0] & (1UL << noABHitInd)){
                                hitPrevFilled = 1;
                                break;
                            }
                        }
                    }
                    if(hitPrevFilled == 0){
                        if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                            //fill spectrum with hit
                            int eGamma = (int)((projAddbackE[projABHitInd]/keVPerBin) + eShift);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[0][eGamma]++;
                                //mcaOut[4][eGamma]++;
                            }
                            //flag that spectrum was filled with this hit
                            for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                                    continue; //skip pileup hit
                                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                                    continue; //skip non-pileup hit
                                }

                                if(projABHitMapping[noABHitInd] == projABHitInd){
                                    specFillFlags[0] |= (1UL << noABHitInd);
                                }
                            }
                        }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                            //time random
                            int eGamma = (int)((projAddbackE[projABHitInd]/keVPerBin) + eShift);
                            if(eGamma>=0 && eGamma<S32K){
                                mcaOut[3][eGamma]++;
                            }
                        }
                    }
                    //check angles of original projection hits, and fill 180 degree spectra
                    //first, check flags
                    uint8_t hitPrevFilledS = 0;
                    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                        if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                            continue; //skip pileup hit
                        }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                            continue; //skip non-pileup hit
                        }

                        if(projABHitMapping[noABHitInd] == projABHitInd){
                            if(specFillFlags[1] & (1UL << noABHitInd)){
                                hitPrevFilledS = 1;
                                break;
                            }
                        }
                    }
                    if(hitPrevFilledS == 0){
                        if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                            for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

                                if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                                    continue; //skip pileup hit
                                }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                                    continue; //skip non-pileup hit
                                }

                                if(!(specFillFlags[1] & (1UL << noABHitInd))){
                                    if(projABHitMapping[noABHitInd] == projABHitInd){
                                        for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

                                            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                                continue; //skip pileup hit
                                            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                                continue; //skip non-pileup hit
                                            }

                                            if(!(specFillFlags[1] & (1UL << noABHitInd2))){
                                                if(((!(projABHitBuildFlags & (1UL << noABHitInd2))) || (projABHitMapping[noABHitInd2] != projABHitInd))&&(!(usedGateHitBuildFlags & (1UL << noABHitInd2)))){
                                                    if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                                        double eGamma = (projAddbackE[projABHitInd]/keVPerBin) + eShift;
                                                        double eGamma2 = 0;
                                                        double tDiff = SUM_TIMING_GATE_MAX + 1000.0; //default value, outside the gate
                                                        if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                            //opposing hit was part of a different projection hit
                                                            eGamma2 = (projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin) + eShift;
                                                            tDiff = (projAddbackT[projABHitMapping[noABHitInd2]] - projAddbackT[projABHitInd]);
                                                        }else{
                                                            //opposing hit was a single non-addback hit
                                                            eGamma2 = (sortedEvt.noABHit[noABHitInd2].energy/keVPerBin) + eShift;
                                                            tDiff = (noABHitTime(&sortedEvt,noABHitInd2) - projAddbackT[projABHitInd]);
                                                        }
                                                        if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                                            //set flags
                                                            for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                if(projABHitMapping[noABHitIndS] == projABHitInd){
                                                                    specFillFlags[1] |= (1UL << noABHitIndS);
                                                                }
                                                            }
                                                            if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                                for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                                    if(projABHitMapping[noABHitIndS] == projABHitMapping[noABHitInd2]){
                                                                        specFillFlags[1] |= (1UL << noABHitIndS);
                                                                    }
                                                                }
                                                            }else{
                                                                specFillFlags[1] |= (1UL << noABHitInd2);
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
                //add in any remaining hits that weren't addback'd
                for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){ //double filling is not in this loop

                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                        continue; //skip pileup hit
                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                        continue; //skip non-pileup hit
                    }

                    if(!(usedGateHitBuildFlags & (1UL << noABHitInd))){
                        double tDiff = (gateAddbackT[gateABHitInd] - noABHitTime(&sortedEvt,noABHitInd));
                        //printf("tDiff: %0.3f\n",tDiff);
                        if(!(specFillFlags[0] & (1UL << noABHitInd))){
                            if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                                //fill spectrum with hit
                                int eGamma = (int)((sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[0][eGamma]++;
                                    //mcaOut[5][eGamma]++;
                                }
                                //flag that spectrum was filled with this hit
                                specFillFlags[0] |= (1UL << noABHitInd);
                                
                            }else if((tDiff >= TRANDOM_GATE_MIN)&&(tDiff <= TRANDOM_GATE_MAX)){
                                //time random
                                int eGamma = (int)((sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift);
                                if(eGamma>=0 && eGamma<S32K){
                                    mcaOut[3][eGamma]++;
                                }
                            }
                        }
                        //check angles of original projection hits, and fill 180 degree spectra
                        if(!(specFillFlags[1] & (1UL << noABHitInd))){
                            if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){ //timing condition
                                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

                                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                                        continue; //skip pileup hit
                                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                                        continue; //skip non-pileup hit
                                    }

                                    if(!(specFillFlags[1] & (1UL << noABHitInd2))){
                                        if(!(usedGateHitBuildFlags & (1UL << noABHitInd2))){
                                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                                double eGamma = (sortedEvt.noABHit[noABHitInd].energy/keVPerBin) + eShift;
                                                double eGamma2 = 0;
                                                double tDiff = SUM_TIMING_GATE_MAX + 1000.0; //default value, outside the gate
                                                if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                    //opposing hit was part of a different projection hit
                                                    eGamma2 = (projAddbackE[projABHitMapping[noABHitInd2]]/keVPerBin) + eShift;
                                                    tDiff = (projAddbackT[projABHitMapping[noABHitInd2]] - noABHitTime(&sortedEvt,noABHitInd));
                                                }else{
                                                    //opposing hit was a single non-addback hit
                                                    eGamma2 = (sortedEvt.noABHit[noABHitInd2].energy/keVPerBin) + eShift;
                                                    tDiff = (noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
                                                }
                                                if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                                    //set flags
                                                    specFillFlags[1] |= (1UL << noABHitInd);
                                                    if(projABHitBuildFlags & (1UL << noABHitInd2)){
                                                        for(int noABHitIndS = 0; noABHitIndS < sortedEvt.header.numNoABHits; noABHitIndS++){
                                                            if(projABHitMapping[noABHitIndS] == projABHitMapping[noABHitInd2]){
                                                                specFillFlags[1] |= (1UL << noABHitIndS);
                                                            }
                                                        }
                                                    }else{
                                                        specFillFlags[1] |= (1UL << noABHitInd2);
                                                    }
                                                    //fill sum spectra
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
                //need to continue to later iterations of the loop, in case there are multiple coincident pairs
            }
            /*if(specFillFlags[0] != 0){
                printf("WARNING: spectrum not filled (gate addback), but fill flags set!\n");
            }*/
        }

        

        if (jentry % 9713 == 0)
            cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;

    } // analysis tree

    cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
    
    fclose(inp);

    return sentries;
  
}

int main(int argc, char **argv){

    const char *sfile;
    const char *outfile;
    double keVPerBin = 1.0;
    double eLow, eHigh, eCorr, gateABRad, projABRad;
    uint8_t forwardPos = 0;
    uint8_t sameCloverOnly = 0;
    uint8_t discardPileup = 0;
    printf("Starting EEGamma_modAB_comptonCorr_mca_SMOL\n");

    if((argc != 10)&&(argc != 11)&&(argc != 12)){
        cout << "Generates dmca gamma-gated Compton background spectra, to be used in background subtraction of Compton scatter features." << endl;
        cout << "Arguments: EEGamma_modAB_comptonCorr_mca_SMOL smol_file EGateLow EGateHigh ECorr gateABRadius projABRadius same_clover_only forward_pos output_dmca_file keV_per_bin discard_pileup" << endl;
        cout << "  *smol_file* can be a single SMOL tree (extension .smol), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
        cout << "  *ECorr* is the energy in keV for which the Compton background is being generated (ie. the original gate energy which the background will be subtracted from)." << endl;
        cout << "  *gateABRadius* and *projABRadius* are in mm. Values of zero correspond to no addback." << endl;
        cout << "  *same_clover_only* should be 1 if the addback will only consider hits in the same clover, 0 otherwise." << endl;
        cout << "  *forward_pos* should be 1 if the detectors are in the forward (110 mm) position, 0 otherwise." << endl;
        cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
        cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
        return 0;
    }else{
        sfile = argv[1];
        eLow = atof(argv[2]);
        eHigh = atof(argv[3]);
        eCorr = atof(argv[4]);
        gateABRad = atof(argv[5]);
        projABRad = atof(argv[6]);
        sameCloverOnly = (uint8_t)(atoi(argv[7]) == 1);
        forwardPos = (uint8_t)(atoi(argv[8]) == 1);
        outfile = argv[9];
        if(argc > 10){
            keVPerBin = atof(argv[10]);
            if(argc > 11){
                discardPileup = (uint8_t)(atoi(argv[11]));
            }
        }
    }

    if(keVPerBin <= 0.0){
        cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
        return 0;
    }

    if(discardPileup > 2){
        cout << "ERROR: Invalid discard_pileup value (" << discardPileup << ")!" << endl;
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

    if(discardPileup == 1){
        printf("Discarding pileup hits.\n");
    }else if(discardPileup == 2){
        printf("Only taking pileup hits.\n");
    }

    const double eShift = (eLow + eHigh)/2.0 - eCorr; //energy to shift all hits by

    uint64_t numSepEvts = 0U;
    if(strcmp(dot + 1, "smol") == 0){
        printf("SMOL tree: %s\nEnergy gate: [%0.2f %0.2f]\nCompton correction energy: %0.2f keV\nGate addback radius: %0.2f mm\nProjection addback radius: %0.2f mm\nDetector position: %s\nClover addback only: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, eCorr, gateABRad, projABRad, (forwardPos == 1) ? "forward" : "back", (sameCloverOnly == 1) ? "yes" : "no", outfile, keVPerBin);
        numSepEvts += SortData(sfile, eLow, eHigh, eShift, gateABRad, projABRad, sameCloverOnly, forwardPos, keVPerBin, discardPileup);
    }else if(strcmp(dot + 1, "list") == 0){
        printf("SMOL tree list: %s\nEnergy gate: [%0.2f %0.2f]\nCompton correction energy: %0.2f keV\nGate addback radius: %0.2f mm\nProjection addback radius: %0.2f mm\nDetector position: %s\nClover addback only: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, eLow, eHigh, eCorr, gateABRad, projABRad, (forwardPos == 1) ? "forward" : "back", (sameCloverOnly == 1) ? "yes" : "no", outfile, keVPerBin);
        
        FILE *listfile;
        char str[256];

        if((listfile=fopen(sfile,"r"))==NULL){
            cout << "ERROR: Cannot open the list file: " << sfile << endl;
            return 0;
        }else{
            while(!(feof(listfile))){//go until the end of file is reached
                if(fgets(str,256,listfile)!=NULL){ //get an entire line
                    str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
                    numSepEvts += SortData(str, eLow, eHigh, eShift, gateABRad, projABRad, sameCloverOnly, forwardPos, keVPerBin, discardPileup);
                }
            }
        }
    }else{
        cout << "ERROR: improper file extension for SMOL tree or list (should be .smol or .list)." << endl;
        return 0;
    }
    
    WriteData(outfile);
    cout << "Wrote " << numSepEvts << " separated events to: " << outfile << endl;

    return 0;
}
