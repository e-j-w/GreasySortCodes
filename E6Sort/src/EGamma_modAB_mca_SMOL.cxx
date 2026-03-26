//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EGamma_modAB_mca_SMOL_cxx
#include "common.cxx"
#include "EGamma_modAB_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

uint8_t ABHitMapping[64]; //arrays specifying which hits correspond to which addback hits
double addbackE[64];
double addbackT[64];

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

uint64_t SortData(const char *sfile, const double projABRad, const uint8_t sameCloverOnly, const uint8_t forwardPos, const double keVPerBin, const uint8_t discardPileup){

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

        //reset flags
        uint64_t abHitBuildFlags = 0;
        uint8_t numProjABHitsBuilt = 0;
        memset(ABHitMapping,0,sizeof(ABHitMapping));
        memset(addbackE,0,sizeof(addbackE));
        memset(addbackT,0,sizeof(addbackT));

        //build addback hits, using the addback radius,
        //and using ABHitMapping to track which hits are grouped together
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
                                if(getGRIFFINHitDistance(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos) < projABRad){
                                    double tDiff = (noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                    if(fabs(tDiff) <= ADDBACK_TIMING_GATE){ //timing condition
                                        if(!(abHitBuildFlags & (1UL << noABHitInd))){
                                            //first hit not yet flagged
                                            if(!(abHitBuildFlags & (1UL << noABHitInd2))){
                                                //neither hit flagged
                                                //build initial addback hit
                                                addbackT[numProjABHitsBuilt] = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                                                if(!(abHitBuildFlags & (1UL << noABHitInd))){
                                                    addbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd].energy;
                                                }
                                                if(!(abHitBuildFlags & (1UL << noABHitInd2))){
                                                    addbackE[numProjABHitsBuilt] += sortedEvt.noABHit[noABHitInd2].energy;
                                                }
                                                //flag hits
                                                abHitBuildFlags |= (1UL << noABHitInd);
                                                abHitBuildFlags |= (1UL << noABHitInd2);
                                                ABHitMapping[noABHitInd] = numProjABHitsBuilt;
                                                ABHitMapping[noABHitInd2] = numProjABHitsBuilt;
                                                abHitBuilt = 1;
                                            }else{
                                                //first hit not flagged, second hit flagged
                                                //add first hit to second hit's addback data
                                                addbackE[ABHitMapping[noABHitInd2]] += sortedEvt.noABHit[noABHitInd].energy;
                                                //flag hit
                                                abHitBuildFlags |= (1UL << noABHitInd);
                                                ABHitMapping[noABHitInd] = ABHitMapping[noABHitInd2];
                                            }
                                        }else{
                                            //first hit flagged
                                            if(!(abHitBuildFlags & (1UL << noABHitInd2))){
                                                //first hit flagged, second not flagged
                                                //add second hit to first hit's addback data
                                                addbackE[ABHitMapping[noABHitInd]] += sortedEvt.noABHit[noABHitInd2].energy;
                                                //flag hit
                                                abHitBuildFlags |= (1UL << noABHitInd2);
                                                ABHitMapping[noABHitInd2] = ABHitMapping[noABHitInd];
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
        //printf("num original hits: %u\n",sortedEvt.header.numNoABHits);
        
        //add hits to the output spectrum
        for(uint8_t projABHitInd=0; projABHitInd < numProjABHitsBuilt; projABHitInd++){
            int eGamma = (int)(addbackE[projABHitInd]/keVPerBin);
            if(eGamma>=0 && eGamma<S32K){
                mcaOut[0][eGamma]++;
            }
        }

        //add remaining non-addbacked hits
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                continue; //skip pileup hit
            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                continue; //skip non-pileup hit
            }

            if(!(abHitBuildFlags & (1UL << noABHitInd))){
                int eGamma = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                if(eGamma>=0 && eGamma<S32K){
                    mcaOut[0][eGamma]++;
                }
            }
        }

        //fill 180 degree sum and projection spectra
        uint64_t usedHits = 0;
        uint64_t usedABHits1 = 0;
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

            if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7))){
                continue; //skip pileup hit
            }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd].core & ((uint8_t)(1) << 7)))){
                continue; //skip non-pileup hit
            }

            if((!(usedHits & (1UL << noABHitInd))) && ((!(abHitBuildFlags & (1UL << noABHitInd))) || (!(usedABHits1 & (1UL << ABHitMapping[noABHitInd]))))){
                int eGamma1 = 0;
                double tGamma1 = 0.0;
                if(abHitBuildFlags & (1UL << noABHitInd)){
                    //hit was part of an addback hit
                    eGamma1 = (int)(addbackE[ABHitMapping[noABHitInd]]/keVPerBin);
                    tGamma1 = addbackT[ABHitMapping[noABHitInd]];
                }else{
                    //hit was a single non-addback hit
                    eGamma1 = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                    tGamma1 = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                }
                uint64_t usedABHits2 = 0;
                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){

                    if((discardPileup == 1) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7))){
                        continue; //skip pileup hit
                    }else if((discardPileup == 2) && (!(sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)(1) << 7)))){
                        continue; //skip non-pileup hit
                    }

                    //if(noABHitInd != noABHitInd2){
                        if((!(usedHits & (1UL << noABHitInd2))) && ((!(abHitBuildFlags & (1UL << noABHitInd2))) || (!(usedABHits2 & (1UL << ABHitMapping[noABHitInd2]))))){
                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                int eGamma2 = 0;
                                double tDiff = SUM_TIMING_GATE_MAX + 1000.0; //default value, outside the gate
                                if(abHitBuildFlags & (1UL << noABHitInd2)){
                                    //opposing hit was part of an addback hit
                                    eGamma2 = (int)(addbackE[ABHitMapping[noABHitInd2]]/keVPerBin);
                                    tDiff = (addbackT[ABHitMapping[noABHitInd2]] - tGamma1);
                                }else{
                                    //opposing hit was a single non-addback hit
                                    eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                    tDiff = (sortedEvt.noABHit[noABHitInd2].timeOffsetNs - tGamma1);
                                }
                                if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                    usedHits |= (1UL << noABHitInd);
                                    usedHits |= (1UL << noABHitInd2);
                                    if(abHitBuildFlags & (1UL << noABHitInd)){
                                        usedABHits1 |= (1UL << ABHitMapping[noABHitInd]);
                                    }
                                    if(abHitBuildFlags & (1UL << noABHitInd2)){
                                        usedABHits2 |= (1UL << ABHitMapping[noABHitInd2]);
                                    }
                                    int eGammaSum = eGamma1 + eGamma2;
                                    if(eGammaSum>=0 && eGammaSum<S32K){
                                        mcaOut[2][eGammaSum]++;
                                    }
                                    if(eGamma1>=0 && eGamma1<S32K){
                                        mcaOut[1][eGamma1]++;
                                    }
                                    if(eGamma2>=0 && eGamma2<S32K){
                                        mcaOut[1][eGamma2]++;
                                    }
                                }
                            }
                        }
                    //}
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

    const char *sfile;
    const char *outfile;
    double keVPerBin = 1.0;
    double projABRad;
    uint8_t forwardPos = 0;
    uint8_t sameCloverOnly = 0;
    uint8_t discardPileup = 0;
    printf("Starting EGamma_modAB_mca_SMOL\n");

    if((argc != 6)&&(argc != 7)&&(argc != 8)){
        cout << "Generates dmca singles spectra." << endl;
        cout << "Arguments: EGamma_modAB_mca_SMOL smol_file projABRadius same_clover_only forward_pos output_dmca_file keV_per_bin discard_pileup" << endl;
        cout << "  *smol_file* can be a single SMOL tree (extension .smol), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
        cout << "  *projABRadius* is in mm. A value of zero corresponds to no addback" << endl;
        cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
        cout << "  *discard_pileup* can be either 0 (false, default if not specified), 1 (true), or 2 (only use pileup hits)." << endl;
        return 0;
    }else{
        sfile = argv[1];
        projABRad = atof(argv[2]);
        sameCloverOnly = (uint8_t)(atoi(argv[3]) == 1);
        forwardPos = (uint8_t)(atoi(argv[4]) == 1);
        outfile = argv[5];
        if(argc > 6){
            keVPerBin = atof(argv[6]);
            if(argc > 7){
                discardPileup = (uint8_t)(atoi(argv[7]));
            }
        }
    }

    if(keVPerBin <= 0.0){
        cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
        return 0;
    }

    if(projABRad < 0.0){
        cout << "ERROR: addback radius must be 0 or larger!" << endl;
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

    uint64_t numSepEvts = 0U;
    if(strcmp(dot + 1, "smol") == 0){
        printf("SMOL tree: %s\nAddback radius: %0.2f mm\nDetector position: %s\nClover addback only: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, projABRad, (forwardPos == 1) ? "forward" : "back", (sameCloverOnly == 1) ? "yes" : "no", outfile, keVPerBin);
        numSepEvts += SortData(sfile, projABRad, sameCloverOnly, forwardPos, keVPerBin, discardPileup);
    }else if(strcmp(dot + 1, "list") == 0){
        printf("SMOL tree list: %s\nAddback radius: %0.2f mm\nDetector position: %s\nClover addback only: %s\nOutput file: %s\n%0.2f keV per bin\n", sfile, projABRad, (forwardPos == 1) ? "forward" : "back", (sameCloverOnly == 1) ? "yes" : "no", outfile, keVPerBin);
        
        FILE *listfile;
        char str[256];

        if((listfile=fopen(sfile,"r"))==NULL){
            cout << "ERROR: Cannot open the list file: " << sfile << endl;
            return 0;
        }else{
            while(!(feof(listfile))){//go until the end of file is reached
                if(fgets(str,256,listfile)!=NULL){ //get an entire line
                    str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
                    numSepEvts += SortData(str, projABRad, sameCloverOnly, forwardPos, keVPerBin, discardPileup);
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
