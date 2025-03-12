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

void EGamma_modAB_mca_SMOL::WriteData(const char* outName){

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

uint64_t EGamma_modAB_mca_SMOL::SortData(const char *sfile, const double projABRad, const double keVPerBin){

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
                uint8_t abHitBuilt = 0;
                for(int noABHitInd2 = 0; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                    if(noABHitInd2 < 64){
                        if(noABHitInd2 != noABHitInd){
                            //check if hits are in neighbouring crystals
                            if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < projABRad){ //FORWARD POSITION (11 cm)
                                double tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                if(tDiff <= ADDBACK_TIMING_GATE){ //timing condition
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
                    //if(noABHitInd != noABHitInd2){
                        if((!(usedHits & (1UL << noABHitInd2))) && ((!(abHitBuildFlags & (1UL << noABHitInd2))) || (!(usedABHits2 & (1UL << ABHitMapping[noABHitInd2]))))){
                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core][sortedEvt.noABHit[noABHitInd2].core] != 0){
                                int eGamma2 = 0;
                                double tDiff = SUM_TIMING_GATE + 1000.0; //default value, outside the gate
                                if(abHitBuildFlags & (1UL << noABHitInd2)){
                                    //opposing hit was part of an addback hit
                                    eGamma2 = (int)(addbackE[ABHitMapping[noABHitInd2]]/keVPerBin);
                                    tDiff = fabs(addbackT[ABHitMapping[noABHitInd2]] - tGamma1);
                                }else{
                                    //opposing hit was a single non-addback hit
                                    eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                    tDiff = fabs(sortedEvt.noABHit[noABHitInd2].timeOffsetNs - tGamma1);
                                }
                                if(tDiff<= SUM_TIMING_GATE){ //timing condition
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

  EGamma_modAB_mca_SMOL *mysort = new EGamma_modAB_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double projABRad;
  printf("Starting EGamma_modAB_mca_SMOL\n");

  if((argc != 4)&&(argc != 5)){
    cout << "Generates dmca singles spectra." << endl;
    cout << "Arguments: EGamma_modAB_mca_SMOL smol_file projABRadius output_dmca_file keV_per_bin" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smole6), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *projABRadius* is in mm. A value of zero corresponds to no addback" << endl;
    cout << "  *keV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    projABRad = atof(argv[2]);
    outfile = argv[3];
    if(argc > 4){
      keVPerBin = atof(argv[4]);
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

  uint64_t numSepEvts = 0U;
  if(strcmp(dot + 1, "smole6") == 0){
    printf("SMOL tree: %s\nAddback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, projABRad, outfile, keVPerBin);
    numSepEvts += mysort->SortData(sfile, projABRad, keVPerBin);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nAddback radius: %0.2f mm\nOutput file: %s\n%0.2f keV per bin\n", sfile, projABRad, outfile, keVPerBin);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, projABRad, keVPerBin);
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
