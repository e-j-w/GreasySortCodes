//use the Makefile!

#define DecayCurveS_cxx
#include "common.cxx"
#include "DecayCurveSMOL.h"

using namespace std;

sorted_evt sortedEvt;

void SortData(const char *sfile, const uint64_t startNumSec)
{
    
  FILE *inp = fopen(sfile, "rb");
  if(inp==NULL){
    printf("ERROR: couldn't open file: %s\n",sfile);
    exit(-1);
  }
  
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
    printf("Total hits:     %Lu\n",totalHits);
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

  printf("Sorting data from file %s, offset by %lu seconds.\n",sfile,startNumSec);
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

    for(int ABpos = 0; ABpos < NGRIFPOS; ABpos++){

      if(addbackT[ABpos] < 0.0){
        //no addback hit in this clover
        continue;
      }

      if(addbackE[ABpos] > MIN_HPGE_EAB){
        tSec = (addbackT[ABpos]/(1.0E9)) + (double)startNumSec;
        //if((jentry % 10000)==0) printf("evtTimeNs: %f, tSec: %f\n",sortedEvt.header.evtTimeNs,tSec);
        hpgeE_time->Fill(addbackE[ABpos], tSec/60.0);
        counts_time->Fill(tSec/60.0);
        if((addbackE[ABpos] >= 676.0)&&(addbackE[ABpos] <= 695.0)){
          counts_time_685->Fill(tSec/60.0);
        }
        if((addbackE[ABpos] >= 930.0)&&(addbackE[ABpos] <= 938.0)){
          counts_time_934->Fill(tSec/60.0);
        }
        if((addbackE[ABpos] >= 1469.0)&&(addbackE[ABpos] <= 1487.0)){
          counts_time_1477->Fill(tSec/60.0);
        }
        counts_time->Fill(tSec/60.0);
      }
    }
    
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        tSec = (noABHitTime(&sortedEvt, noABHitInd)/(1.0E9)) + (double)startNumSec;
        //if((jentry % 10000)==0) printf("evtTimeNs: %f, tSec: %f\n",sortedEvt.header.evtTimeNs,tSec);
        hpgeE_time->Fill(sortedEvt.noABHit[noABHitInd].energy, tSec/60.0);
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

  if(argc == 1){
    cout << "Code sorts decay curve data." << endl;
    cout << "Arguments: DecayCurveSMOL list_file output_file" << endl;
    cout << "The list file should contain two columns (space-delimited) with the SMOL tree filenames and the times (in seconds) at the start of each run." << endl;
    return 0;
  }else if(argc == 3){
    lsfile = argv[1];
    outfile = argv[2];
    printf("List file: %s\nOutput file: %s\n", lsfile, outfile);
  }else{
    printf("ERROR: too many arguments!\nArguments: DecayCurve SMOL smol_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  const char *dot = strrchr(lsfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get list file name." << endl;
    return 0;
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
                        SortData(fileName,(uint64_t)strtol(tok,NULL,10));
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

  myfile->Write();
  myfile->Close();

  return 0;
}
