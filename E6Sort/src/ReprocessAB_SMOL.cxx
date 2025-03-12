//Sort code to re-do addback in SMOL files

#define ReprocessAB_SMOL_cxx
#include "common.h"
#include "evt_fmt.h"
#include "ReprocessAB_SMOL.h"

using namespace std;

FILE *inp, *out;
uint8_t ABHitMapping[64]; //array specifying which hits correspond to which addback hits

uint64_t ReprocessAB_SMOL::SortData(){

  uint64_t sentries = 0U;
  sorted_evt sortedEvt;
  uint8_t footerVal = 227U;
  uint64_t hitBuildFlags = 0;
  uint64_t numSeparatedEvents = 0;

  fread(&sentries,sizeof(uint64_t),1,inp);

  printf("\nReprocessing events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event from input file
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    //reset flags
    hitBuildFlags = 0;
    memset(ABHitMapping,0,sizeof(ABHitMapping));

    //reset addback in event
    sortedEvt.header.numABHits = 0;
    memset(sortedEvt.ABHit,0,sizeof(sortedEvt.ABHit));

    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(noABHitInd < 64){
        uint8_t abHitBuilt = 0;
        for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
          if(noABHitInd2 < 64){
            //check if hits are in neighbouring crystals
            //if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < 150.0){ //FORWARD POSITION (11 cm)
            //if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < 80.0){ //FORWARD POSITION (11 cm)
            if(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1) < 65.0){ //FORWARD POSITION (11 cm)
              float tDiff = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2);
              if((tDiff >= hpgehpgeABGate[0])&&(tDiff <= hpgehpgeABGate[1])){ //timing condition
                if(!(hitBuildFlags & (1UL << noABHitInd))){
                  //first hit not yet flagged
                  if(!(hitBuildFlags & (1UL << noABHitInd2))){
                    //neither hit flagged
                    //build initial addback hit
                    sortedEvt.ABHit[sortedEvt.header.numABHits].core = sortedEvt.noABHit[noABHitInd].core;
                    sortedEvt.ABHit[sortedEvt.header.numABHits].timeOffsetNs = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
                    if(!(hitBuildFlags & (1UL << noABHitInd))){
                      sortedEvt.ABHit[sortedEvt.header.numABHits].energy += sortedEvt.noABHit[noABHitInd].energy;
                    }
                    if(!(hitBuildFlags & (1UL << noABHitInd2))){
                      sortedEvt.ABHit[sortedEvt.header.numABHits].energy += sortedEvt.noABHit[noABHitInd2].energy;
                    }
                    //flag hits
                    hitBuildFlags |= (1UL << noABHitInd);
                    hitBuildFlags |= (1UL << noABHitInd2);
                    ABHitMapping[noABHitInd] = sortedEvt.header.numABHits;
                    ABHitMapping[noABHitInd2] = sortedEvt.header.numABHits;
                    abHitBuilt = 1;
                  }else{
                    //first hit not flagged, second hit flagged
                    //add first hit to second hit's addback data
                    sortedEvt.ABHit[ABHitMapping[noABHitInd2]].energy += sortedEvt.noABHit[noABHitInd].energy;
                    //flag hit
                    hitBuildFlags |= (1UL << noABHitInd);
                    ABHitMapping[noABHitInd] = ABHitMapping[noABHitInd2];
                  }
                }else{
                  //first hit flagged
                  if(!(hitBuildFlags & (1UL << noABHitInd2))){
                    //first hit flagged, second not flagged
                    //add second hit to first hit's addback data
                    sortedEvt.ABHit[ABHitMapping[noABHitInd]].energy += sortedEvt.noABHit[noABHitInd2].energy;
                    //flag hit
                    hitBuildFlags |= (1UL << noABHitInd2);
                    ABHitMapping[noABHitInd2] = ABHitMapping[noABHitInd];
                  }
                }
              }
            }
          }
        }
        if(abHitBuilt != 0){
          //addback hit was just built
          sortedEvt.header.numABHits++;
          if((sortedEvt.header.numABHits>=MAXNUMHPGEHIT)||(sortedEvt.header.numABHits>=MAX_EVT_HIT)){
            break;
          }
        }
      }
      
    }

    //copy over any non-addbacked hits as well
    if(!((sortedEvt.header.numABHits>=MAXNUMHPGEHIT)||(sortedEvt.header.numABHits>=MAX_EVT_HIT))){
      for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
        if(!(hitBuildFlags & (1UL << noABHitInd))){ //hit wasn't flagged for addback
          sortedEvt.ABHit[sortedEvt.header.numABHits].timeOffsetNs = sortedEvt.noABHit[noABHitInd].timeOffsetNs;
          sortedEvt.ABHit[sortedEvt.header.numABHits].core = sortedEvt.noABHit[noABHitInd].core;
          sortedEvt.ABHit[sortedEvt.header.numABHits].energy = sortedEvt.noABHit[noABHitInd].energy;
          sortedEvt.header.numABHits++;
          if((sortedEvt.header.numABHits>=MAXNUMHPGEHIT)||(sortedEvt.header.numABHits>=MAX_EVT_HIT)){
            break;
          }
        }
      }
    }
    

    //write out data
    fwrite(&sortedEvt.header,sizeof(evt_header),1,out);
    //write hits, without segment data (no segments for GRIFFIN)
    for(int i = 0; i<sortedEvt.header.numABHits;i++){
      fwrite(&sortedEvt.ABHit[i].timeOffsetNs,sizeof(float),1,out);
      fwrite(&sortedEvt.ABHit[i].energy,sizeof(float),1,out);
      fwrite(&sortedEvt.ABHit[i].core,sizeof(uint8_t),1,out);
    }
    for(int i = 0; i<sortedEvt.header.numNoABHits;i++){
      fwrite(&sortedEvt.noABHit[i].timeOffsetNs,sizeof(float),1,out);
      fwrite(&sortedEvt.noABHit[i].energy,sizeof(float),1,out);
      fwrite(&sortedEvt.noABHit[i].core,sizeof(uint8_t),1,out);
    }
    //write footer value
    fwrite(&footerVal,sizeof(uint8_t),1,out);

    numSeparatedEvents++;

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Addback reprocessing complete" << endl;
  cout << "Number of separated events: " << numSeparatedEvents << " ("<< 100.0f*numSeparatedEvents/(float)sentries << "\% of total)" << endl;

  return (uint64_t)numSeparatedEvents;

}
int main(int argc, char **argv){

  ReprocessAB_SMOL *mysort = new ReprocessAB_SMOL();

  char const *sfile;
  char const *soutfile;

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Arguments: ReprocessAB_SMOL smol_file output_smolfile" << endl;
    cout << "*smol_file* can be a single file (extension .smol)." << endl;
    cout << "Default values will be used if arguments (other than analysis_tree) are omitted." << endl;
    return 0;
  }else if (argc == 2){
    sfile = argv[1];
    soutfile = "Output.smole6";
  }else if (argc == 3){
    sfile = argv[1];
    soutfile = argv[2];
  }else{
    printf("Incorrect arguments\nArguments: ReprocessAB_SMOL smol_file output_smolfile\n");
    return 0;
  }

  printf("Starting ReprocessAB_SMOL code\n");
  if(strcmp(sfile,soutfile)==0){
    cout << "ERROR: input and output filenames are the same!" << endl;
    return 0;
  }

  cout << "Input SMOL file: " << sfile << endl;
  cout << "Output SMOL file: " << soutfile << endl;

  //setup the input file
  inp = fopen(sfile, "r");
  if(inp == NULL){
    cout << "ERROR: couldn't open file " << sfile << endl;
    return 0;
  }

  //setup the output file
  out = fopen(soutfile, "wb");
  if(out == NULL){
    cout << "ERROR: couldn't open output file " << soutfile << endl;
    return 0;
  }

  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);
  numSepEvts += mysort->SortData();

  //write the number of separated events to the beginning of the file
  fseek(out,0,SEEK_SET);
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);
  cout << "Wrote " << numSepEvts << " separated events to: " << soutfile << endl;
  fclose(out);

  return 0;
}
