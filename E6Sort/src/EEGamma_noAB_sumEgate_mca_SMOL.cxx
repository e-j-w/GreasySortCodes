//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEGamma_noAB_sumEgate_mca_SMOL_cxx
#include "common.cxx"
#include "EEGamma_noAB_sumEgate_mca_SMOL.h"

using namespace std;

uint8_t hitMap180deg[64][64]; //1st index = crystal of hit, 2nd index = crystal of 2nd hit, val = 1 indicates 180 degree summing occurs

void EEGamma_noAB_sumEgate_mca_SMOL::WriteData(const char* outName){

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

uint64_t EEGamma_noAB_sumEgate_mca_SMOL::SortData(const char *sfile, const double eLow, const double eHigh, const double keVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

    //construct 180 degree summing hit map
    memset(hitMap180deg,0,sizeof(hitMap180deg));
    for(uint8_t i=0;i<64;i++){ //first core
        for(uint8_t j=0;j<64;j++){ //coinc core
            if(i!=j){
                if(getGRIFFINVector(i,1).Angle(getGRIFFINVector(j,1))*180.0/PI > 175.0){
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

        //There was a double counting issue that could only be fixed by dropping hits with identical energy in the same event.
        //I'm not sure why this works...
        //double counting window width is affected by size of gate (eHigh-eLow)
        //double counting was not fixed by:
        //-restricting number of hits to 2
        //-starting noABHitInd2 at 0
        //-requiring first hit to be higher energy than the 2nd
        //-changing # of keV per bin
        //-disabling filling of 180 degree sum/projection spectra
        //-breaking out of noABHitInd2 loop after filling
        //-disabling filling of 1st gamma energy
        //-disabling filling of 2nd gamma energy

        uint64_t hitFlags = 0;
        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
            if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
                for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
                    if((!(hitFlags & (1U << noABHitInd))) && (!(hitFlags & (1U << noABHitInd2)))){
                        if(sortedEvt.noABHit[noABHitInd2].energy > MIN_HPGE_EAB){
                            double sumE = (double)(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy);
                            //printf("sumE: %0.3f\n",sumE);
                            if((sumE >= eLow)&&(sumE <= eHigh)){
                                Double_t tDiff = fabs(noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2));
                                if(tDiff <= ADDBACK_TIMING_GATE){
                                    int eGamma1 = (int)(sortedEvt.noABHit[noABHitInd].energy/keVPerBin);
                                    int eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                    if(eGamma1>=0 && eGamma1<S32K){
                                        mcaOut[0][eGamma1]++;
                                    }
                                    if(eGamma2>=0 && eGamma2<S32K){
                                        mcaOut[0][eGamma2]++;
                                    }
                                    hitFlags |= (1U << noABHitInd);
                                    hitFlags |= (1U << noABHitInd2);

                                    //printf("Entry %li, filled at inds: %i %i\n",jentry,noABHitInd,noABHitInd2);
                                    /*if(fabs( sortedEvt.noABHit[noABHitInd].energy - sortedEvt.noABHit[noABHitInd2].energy) < 15.0f){
                                        printf("\ninds: %i %i, tDiff: %0.3f\n",noABHitInd,noABHitInd2,tDiff);
                                        for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                                            printf(" Gamma %i, E: %0.3f, core: %u, hit flag: %u\n",noABHitInd3,sortedEvt.noABHit[noABHitInd3].energy,sortedEvt.noABHit[noABHitInd3].core & 63U,hitFlags & (1U << noABHitInd3));
                                        }
                                    }*/
                                    for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                                        if((noABHitInd3 != noABHitInd)&&((noABHitInd3 != noABHitInd2))){
                                            if(hitMap180deg[sortedEvt.noABHit[noABHitInd3].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] != 0){
                                                tDiff = (noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd2));
                                                if((tDiff >= SUM_TIMING_GATE_MIN)&&(tDiff <= SUM_TIMING_GATE_MAX)){ //timing condition
                                                    int eGamma3 = (int)(sortedEvt.noABHit[noABHitInd3].energy/keVPerBin);
                                                    int eGammaSum = (int)((sortedEvt.noABHit[noABHitInd3].energy + sortedEvt.noABHit[noABHitInd2].energy)/keVPerBin);
                                                    if(eGammaSum>=0 && eGammaSum<S32K){
                                                        mcaOut[2][eGammaSum]++;
                                                    }
                                                    if(eGamma3>=0 && eGamma3<S32K){
                                                        mcaOut[1][eGamma3]++;
                                                    }
                                                    int eGamma2 = (int)(sortedEvt.noABHit[noABHitInd2].energy/keVPerBin);
                                                    if(eGamma2>=0 && eGamma2<S32K){
                                                        mcaOut[1][eGamma2]++;
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

  EEGamma_noAB_sumEgate_mca_SMOL *mysort = new EEGamma_noAB_sumEgate_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 1.0;
  double eLow, eHigh;
  printf("Starting EEGamma_noAB_sumEgate_mca_SMOL\n");

  if((argc != 5)&&(argc != 6)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEGamma_noAB_sumEgate_mca_SMOL smol_file EGateLow EGateHigh output_dmca_file keV_per_bin" << endl;
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
