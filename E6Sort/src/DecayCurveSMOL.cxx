//use the Makefile!

#define DecayCurveS_cxx
#include "common.cxx"
#include "DecayCurveSMOL.h"

using namespace std;

sorted_evt sortedEvt;

const Double_t gate511[2] = {506.0, 515.0};
const Double_t gate934[2] = {924.0, 944.0};
const Double_t gate955[2] = {945.0, 965.0};
const Double_t gate1238[2] = {1235.0, 1240.5};

void SortData(const char *sfile, const uint64_t startNumSec, const uint8_t forwardPos)
{
    
  FILE *inp = fopen(sfile, "rb");
  if(inp==NULL){
    printf("ERROR: couldn't open file: %s\n",sfile);
    exit(-1);
  }

  printf("Sorting data from file %s, offset by %lu seconds.\n",sfile,startNumSec);
  
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
    printf("Total hits in file:     %Lu\n",totalHits);
    long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
    printf("Fraction of hits with type 1 (no pileup): %Lf\n",frac);
  }
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  sorted_evt sortedEvt;

  /*cout << "TIGRESS positions: " << endl;
  for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
      TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
      cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
  }*/
  
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      printf("ERROR: bad event data in entry %lu.\n",jentry);
      exit(-1);
    }

    //construct addback energies and times
    memset(addbackE,0,sizeof(addbackE));
    memset(maxABHitE,0,sizeof(maxABHitE));
    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){
      addbackT[ABpos] = -1.0; //default value
    }
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

      int ABHitPos = (sortedEvt.noABHit[noABHitInd].core & 63U)/4;

      if(ABHitPos < NGRIFPOS){

        //check timing criteria
        if(addbackT[ABHitPos] >= 0.0){
          //there are hits in this clover
          if(fabs(noABHitTime(&sortedEvt,noABHitInd) - (double)addbackT[ABHitPos]) > ADDBACK_TIMING_GATE){
            //hit not in time coincidence with other hits
            if(sortedEvt.noABHit[noABHitInd].energy > maxABHitE[ABHitPos]){
              //higher energy, not time coincident
              //make this hit the new hit
              addbackT[ABHitPos] = noABHitTime(&sortedEvt,noABHitInd);
              maxABHitE[ABHitPos] = sortedEvt.noABHit[noABHitInd].energy;
              addbackE[ABHitPos] = sortedEvt.noABHit[noABHitInd].energy;
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
        if(sortedEvt.noABHit[noABHitInd].energy > maxABHitE[ABHitPos]){
          addbackT[ABHitPos] = noABHitTime(&sortedEvt,noABHitInd);
          maxABHitE[ABHitPos] = sortedEvt.noABHit[noABHitInd].energy;
        }
        addbackE[ABHitPos] += sortedEvt.noABHit[noABHitInd].energy; 
      }
    }
    
    Double_t tSec = 0.0;
    uint8_t pos1238 = NGRIFPOS; //position index for 1238 keV hit
    uint8_t pos934 = NGRIFPOS; //position index for 934 keV hit
    uint8_t pos955 = NGRIFPOS; //position index for 955 keV hit

    //addback
    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){
      if((pos1238 == NGRIFPOS)&&(addbackE[ABpos] >= gate1238[0])&&(addbackE[ABpos] <= gate1238[1])){
        pos1238 = ABpos;
      }\
      if((pos955 == NGRIFPOS)&&(addbackE[ABpos] >= gate955[0])&&(addbackE[ABpos] <= gate955[1])){
        pos955 = ABpos;
      }
      if((pos934 == NGRIFPOS)&&(addbackE[ABpos] >= gate934[0])&&(addbackE[ABpos] <= gate934[1])){
        pos934 = ABpos;
      }
    }
    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){

      if(addbackT[ABpos] < 0.0){
        //no addback hit in this clover
        continue;
      }

      if(addbackE[ABpos] > MIN_HPGE_EAB){
        tSec = (addbackT[ABpos]/(1.0E9)) + (double)startNumSec;
        //if((jentry % 10000)==0) printf("evtTimeNs: %f, tSec: %f\n",sortedEvt.header.evtTimeNs,tSec);
        addE_time->Fill(addbackE[ABpos], tSec/60.0);
        addE->Fill(addbackE[ABpos]);
        counts_time->Fill(tSec/60.0);
        if((pos934 < NGRIFPOS)&&(pos934 != ABpos)){
          Double_t tDiff934 = addbackT[ABpos] - addbackT[pos934];
          if((tDiff934 >= COINC_TIMING_GATE_MIN)&&(tDiff934 <= COINC_TIMING_GATE_MAX)){
            addE934->Fill(addbackE[ABpos]);
          }
        }
        if((addbackE[ABpos] >= 676.0)&&(addbackE[ABpos] <= 695.0)){
          counts_time_685->Fill(tSec/60.0);
        }
        if((addbackE[ABpos] >= gate934[0])&&(addbackE[ABpos] <= gate934[1])){
          counts_time_934->Fill(tSec/60.0);
        }
        if((addbackE[ABpos] >= 964.0)&&(addbackE[ABpos] <= 984.0)){
          counts_time_934_bg->Fill(tSec/60.0);
        }
        if((addbackE[ABpos] >= 1469.0)&&(addbackE[ABpos] <= 1487.0)){
          counts_time_1477->Fill(tSec/60.0);
        }
        counts_time->Fill(tSec/60.0);
      }
      if((addbackE[ABpos] >= gate511[0])&&(addbackE[ABpos] <= gate511[1])){
        for(int ABpos2 = ABpos+1; ABpos2 < NGRIFPOS; ABpos2++){

          if(addbackT[ABpos2] < 0.0){
            //no addback hit in this clover
            continue;
          }

          Double_t tDiff = addbackT[ABpos2] - addbackT[ABpos];
          if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){
            if(hitMap180degAB[ABpos][ABpos2] == 1){
              addE_time_betaplus->Fill(addbackE[ABpos2], tSec/60.0);
              if((pos934 < NGRIFPOS)&&(pos934 != ABpos)&&(pos934 != ABpos2)){
                Double_t tDiff934 = addbackT[ABpos] - addbackT[pos934];
                Double_t tDiff934_2 = addbackT[ABpos2] - addbackT[pos934];
                if(((tDiff934 >= COINC_TIMING_GATE_MIN)&&(tDiff934 <= COINC_TIMING_GATE_MAX)) || ((tDiff934_2 >= COINC_TIMING_GATE_MIN)&&(tDiff934_2 <= COINC_TIMING_GATE_MAX))){
                  addE_time_betaplus_934->Fill(addbackE[ABpos2], tSec/60.0);
                  if((addbackE[ABpos2] >= gate511[0])&&(addbackE[ABpos2] <= gate511[1])){
                    counts_time_511betaplusAB_934->Fill(tSec/60.0);
                  }
                }
              }
              if((pos955 < NGRIFPOS)&&(pos955 != ABpos)&&(pos955 != ABpos2)){
                Double_t tDiff955 = addbackT[ABpos] - addbackT[pos955];
                Double_t tDiff955_2 = addbackT[ABpos2] - addbackT[pos955];
                if(((tDiff955 >= COINC_TIMING_GATE_MIN)&&(tDiff955 <= COINC_TIMING_GATE_MAX)) || ((tDiff955_2 >= COINC_TIMING_GATE_MIN)&&(tDiff955_2 <= COINC_TIMING_GATE_MAX))){
                  addE_time_betaplus_955->Fill(addbackE[ABpos2], tSec/60.0);
                  if((addbackE[ABpos2] >= gate511[0])&&(addbackE[ABpos2] <= gate511[1])){
                    counts_time_511betaplusAB_955->Fill(tSec/60.0);
                  }
                }
              }
              if((addbackE[ABpos2] >= gate511[0])&&(addbackE[ABpos2] <= gate511[1])){
                counts_time_511betaplusAB->Fill(tSec/60.0);
                if((pos1238 < NGRIFPOS)&&(pos1238 != ABpos)&&(pos1238 != ABpos2)){
                  Double_t tDiff1238 = addbackT[ABpos] - addbackT[pos1238];
                  Double_t tDiff1238_2 = addbackT[ABpos2] - addbackT[pos1238];
                  if(((tDiff1238 >= COINC_TIMING_GATE_MIN)&&(tDiff1238 <= COINC_TIMING_GATE_MAX)) || ((tDiff1238_2 >= COINC_TIMING_GATE_MIN)&&(tDiff1238_2 <= COINC_TIMING_GATE_MAX))){
                    counts_time_511betaplusAB_1238->Fill(tSec/60.0);
                  }
                }
                for(int ABpos3 = 0; ABpos3 < NGRIFPOS; ABpos3++){

                  if(addbackT[ABpos3] < 0.0){
                    //no addback hit in this clover
                    continue;
                  }

                  if((ABpos3 != ABpos)&&(ABpos3 != ABpos2)){
                    Double_t tDiff2 = addbackT[ABpos3] - addbackT[ABpos];
                    Double_t tDiff3 = addbackT[ABpos3] - addbackT[ABpos2];
                    if(((tDiff2 >= COINC_TIMING_GATE_MIN)&&(tDiff2 <= COINC_TIMING_GATE_MAX)) || ((tDiff3 >= COINC_TIMING_GATE_MIN)&&(tDiff3 <= COINC_TIMING_GATE_MAX))){
                      if(addbackE[ABpos3] > MIN_HPGE_EAB){
                        addE_time_betaplus_3rdgamma->Fill(addbackE[ABpos3], tSec/60.0);
                      }
                    }
                  }
                }
              }
            }
            //all 511 coincidences (not just 180 degrees)
            if((addbackE[ABpos2] >= gate511[0])&&(addbackE[ABpos2] <= gate511[1])){
              angle511betaplusAB_time->Fill(getGRIFFINVector(ABpos*4,forwardPos).Angle(getGRIFFINVector(ABpos2*4,forwardPos))*180.0/PI,tSec/60.0);
              for(int ABpos3 = 0; ABpos3 < NGRIFPOS; ABpos3++){

                if(addbackT[ABpos3] < 0.0){
                  //no addback hit in this clover
                  continue;
                }

                if((ABpos3 != ABpos)&&(ABpos3 != ABpos2)){
                  Double_t tDiff2 = addbackT[ABpos3] - addbackT[ABpos];
                  Double_t tDiff3 = addbackT[ABpos3] - addbackT[ABpos2];
                  if(((tDiff2 >= COINC_TIMING_GATE_MIN)&&(tDiff2 <= COINC_TIMING_GATE_MAX)) || ((tDiff3 >= COINC_TIMING_GATE_MIN)&&(tDiff3 <= COINC_TIMING_GATE_MAX))){
                    if(addbackE[ABpos3] > MIN_HPGE_EAB){
                      angle511betaplusAB_ecoinc->Fill(getGRIFFINVector(ABpos*4,forwardPos).Angle(getGRIFFINVector(ABpos2*4,forwardPos))*180.0/PI,addbackE[ABpos3]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    pos1238 = sortedEvt.header.numNoABHits; //reset position index for 1238 keV hit
    pos934 = sortedEvt.header.numNoABHits; //reset position index for 934 keV hit
    
    //non-addback (single crystal)
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if((pos1238 == sortedEvt.header.numNoABHits)&&(sortedEvt.noABHit[noABHitInd].energy >= gate1238[0])&&(sortedEvt.noABHit[noABHitInd].energy <= gate1238[1])){
        pos1238 = noABHitInd;
      }
      if((pos934 == sortedEvt.header.numNoABHits)&&(sortedEvt.noABHit[noABHitInd].energy >= gate934[0])&&(sortedEvt.noABHit[noABHitInd].energy <= gate934[1])){
        pos934 = noABHitInd;
      }
    }
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        tSec = (noABHitTime(&sortedEvt, noABHitInd)/(1.0E9)) + (double)startNumSec;
        //if((jentry % 10000)==0) printf("evtTimeNs: %f, tSec: %f\n",sortedEvt.header.evtTimeNs,tSec);
        hpgeE_time->Fill(sortedEvt.noABHit[noABHitInd].energy, tSec/60.0);
        hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy);
        if((pos934 < sortedEvt.header.numNoABHits)&&(pos934 != noABHitInd)){
          Double_t tDiff934 = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,pos934);
          if((tDiff934 >= COINC_TIMING_GATE_MIN)&&(tDiff934 <= COINC_TIMING_GATE_MAX)){
            hpgeE934->Fill(sortedEvt.noABHit[noABHitInd].energy);
          }
        }
        if((sortedEvt.noABHit[noABHitInd].energy >= gate511[0])&&(sortedEvt.noABHit[noABHitInd].energy <= gate511[1])){
          for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            Double_t tDiff = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd);
            if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){
              if(hitMap180deg[sortedEvt.noABHit[noABHitInd].core & 63U][sortedEvt.noABHit[noABHitInd2].core & 63U] == 1){
                hpgeE_time_betaplus->Fill(sortedEvt.noABHit[noABHitInd2].energy, tSec/60.0);
                if((pos934 < sortedEvt.header.numNoABHits)&&(pos934 != noABHitInd)&&(pos934 != noABHitInd2)){
                  Double_t tDiff934_2 = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,pos934);
                  Double_t tDiff934 = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,pos934);
                  if(((tDiff934 >= COINC_TIMING_GATE_MIN)&&(tDiff934 <= COINC_TIMING_GATE_MAX)) || ((tDiff934_2 >= COINC_TIMING_GATE_MIN)&&(tDiff934_2 <= COINC_TIMING_GATE_MAX))){
                    hpgeE_time_betaplus_934->Fill(sortedEvt.noABHit[noABHitInd2].energy, tSec/60.0);
                    if((sortedEvt.noABHit[noABHitInd2].energy >= gate511[0])&&(sortedEvt.noABHit[noABHitInd2].energy <= gate511[1])){
                      counts_time_511betaplus_934->Fill(tSec/60.0);
                    }
                  }
                }
                
                if((sortedEvt.noABHit[noABHitInd2].energy >= gate511[0])&&(sortedEvt.noABHit[noABHitInd2].energy <= gate511[1])){
                  counts_time_511betaplus->Fill(tSec/60.0);
                  if((pos1238 < sortedEvt.header.numNoABHits)&&(pos1238 != noABHitInd)&&(pos1238 != noABHitInd2)){
                    Double_t tDiff1238_2 = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,pos1238);
                    Double_t tDiff1238 = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,pos1238);
                    if(((tDiff1238 >= COINC_TIMING_GATE_MIN)&&(tDiff1238 <= COINC_TIMING_GATE_MAX)) || ((tDiff1238_2 >= COINC_TIMING_GATE_MIN)&&(tDiff1238_2 <= COINC_TIMING_GATE_MAX))){
                      counts_time_511betaplus_1238->Fill(tSec/60.0);
                    }
                  }
                  for(int noABHitInd3 = 0; noABHitInd3 < NGRIFPOS; noABHitInd3++){
                    if((noABHitInd3 != noABHitInd)&&(noABHitInd3 != noABHitInd2)){
                      Double_t tDiff2 = noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd);
                      Double_t tDiff3 = noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd2);
                      if(((tDiff2 >= COINC_TIMING_GATE_MIN)&&(tDiff2 <= COINC_TIMING_GATE_MAX)) || ((tDiff3 >= COINC_TIMING_GATE_MIN)&&(tDiff3 <= COINC_TIMING_GATE_MAX))){
                        if(sortedEvt.noABHit[noABHitInd3].energy > MIN_HPGE_EAB){
                          addE_time_betaplus_3rdgamma->Fill(sortedEvt.noABHit[noABHitInd3].energy, tSec/60.0);
                        }
                      }
                    }
                  }
                }
              }
              //all 511 coincidences (not just 180 degrees)
              if((sortedEvt.noABHit[noABHitInd2].energy >= gate511[0])&&(sortedEvt.noABHit[noABHitInd2].energy <= gate511[1])){
                angle511betaplus_time->Fill(getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos))*180.0/PI,tSec/60.0);
                for(int noABHitInd3 = 0; noABHitInd3 < NGRIFPOS; noABHitInd3++){
                  if((noABHitInd3 != noABHitInd)&&(noABHitInd3 != noABHitInd2)){
                    Double_t tDiff2 = noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd);
                    Double_t tDiff3 = noABHitTime(&sortedEvt,noABHitInd3) - noABHitTime(&sortedEvt,noABHitInd2);
                    if(((tDiff2 >= COINC_TIMING_GATE_MIN)&&(tDiff2 <= COINC_TIMING_GATE_MAX)) || ((tDiff3 >= COINC_TIMING_GATE_MIN)&&(tDiff3 <= COINC_TIMING_GATE_MAX))){
                      if(sortedEvt.noABHit[noABHitInd3].energy > MIN_HPGE_EAB){
                        angle511betaplus_ecoinc->Fill(getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos))*180.0/PI,sortedEvt.noABHit[noABHitInd3].energy);
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

    if (jentry % 90713 == 0){
      //cout << "tNs: " << sortedEvt.header.evtTimeNs << ", tSec: " << tSec << endl;
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
    }
      
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete" << endl;
  
  fclose(inp);
}
int main(int argc, char **argv)
{

  const char *lsfile;
  const char *outfile;
  printf("Starting DecayCurveSMOL\n");
  uint8_t forwardPos = 0;

  if(argc == 1){
    cout << "Code sorts decay curve data." << endl;
    cout << "Arguments: DecayCurveSMOL list_file output_file forward_pos" << endl;
    cout << "The list file should contain two columns (space-delimited) with the SMOL tree filenames and the times (in seconds) at the start of each run." << endl;
    cout << "*forward_pos* specifies whether the array is in the forward (110 mm) position (optional, default 0=no)." << endl;
    return 0;
  }else if(argc >= 3){
    lsfile = argv[1];
    outfile = argv[2];
    printf("List file: %s\nOutput file: %s\n", lsfile, outfile);
    if(argc > 3){
      forwardPos = (uint8_t)atoi(argv[3]);
    }
  }else{
    printf("ERROR: too many arguments!\nArguments: DecayCurve SMOL smol_file output_file\n");
    return 0;
  }

  if(forwardPos == 0){
    printf("Forward position: no\n");
  }else{
    printf("Forward position: yes\n");
  }

  theApp=new TApplication("App", &argc, argv);

  const char *dot = strrchr(lsfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get list file name." << endl;
    return 0;
  }

  //construct 180 degree hit maps (for 511s from beta+)
  //both for addback and non-addback
  memset(hitMap180deg,0,sizeof(hitMap180deg));
  for(uint8_t i=0;i<64;i++){ //first core
    for(uint8_t j=0;j<64;j++){ //coinc core
        if(i!=j){
            if(getGRIFFINVector(i,forwardPos).Angle(getGRIFFINVector(j,forwardPos))*180.0/PI > 175.0){ //same effect for any value down to 165 degrees
                hitMap180deg[i][j] = 1;
                continue; //check the next coinc core
            }
        }
    }
  }
  memset(hitMap180degAB,0,sizeof(hitMap180degAB));
  for(uint8_t i=0;i<16;i++){ //first core
    for(uint8_t j=0;j<16;j++){ //coinc core
      if(i!=j){
        for(uint8_t i_cry=0; i_cry<4; i_cry++){
          for(uint8_t j_cry=0; j_cry<4; j_cry++){
            if(getGRIFFINVector((i*4)+i_cry,forwardPos).Angle(getGRIFFINVector((j*4)+j_cry,forwardPos))*180.0/PI > 175.0){ //same effect for any value down to 165 degrees
              hitMap180degAB[i][j] = 1;
              //printf("Pair %u and %u are at 180 degrees.\n",i,j);
            }
          }
        }
      }
    }
  }

  if(strcmp(dot + 1, "list") == 0){
    
    FILE *listfile;
    char str[256];
    char fileName[128];
    char *tok;

    if((listfile=fopen(lsfile,"r"))==NULL){
        cout << "ERROR: Cannot open the list file: " << lsfile << endl;
        return 0;
    }else{
        InitialiseHists();
        while(!(feof(listfile))){//go until the end of file is reached
            if(fgets(str,256,listfile)!=NULL){ //get an entire line
                str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
                tok = strtok(str,"\t ");
                if(tok!=NULL){
                    strncpy(fileName,tok,127);
                    tok = strtok(NULL,"");
                    if(tok!=NULL){
                        SortData(fileName,(uint64_t)strtol(tok,NULL,10),forwardPos);
                    }else{
                        printf("ERROR: couldn't get time for run in file: %s\n",fileName);
                        exit(-1);
                    }
                }
            }
        }
    }
  }else{
    cout << "ERROR: improper file extension for input list (should be .list)." << endl;
    return 0;
  }

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *dcdir = myfile->mkdir("DecayCurve");
  dcdir->cd();
  dcList->Write();
  myfile->cd();
  TDirectory *dcABdir = myfile->mkdir("AddbackDecayCurve");
  dcABdir->cd();
  dcABList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();

  return 0;
}
