//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define ETIP_mca_SMOL_cxx
#include "common.h"
#include "ETIP_mca_SMOL.h"

using namespace std;

void ETIP_mca_SMOL::WriteData(const char* outName){

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

void ETIP_mca_SMOL::SortData(char const *sfile, double MeVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    memset(&sortedEvt,0,sizeof(sorted_evt));
    footerVal = 0;
    fread(&sortedEvt.header,sizeof(evt_header),1,inp);
    for(int i = 0; i<sortedEvt.header.numTigHits;i++){
      fread(&sortedEvt.tigHit[i],sizeof(tig_hit),1,inp);
    }
    for(int i = 0; i<sortedEvt.header.numCsIHits;i++){
      fread(&sortedEvt.csiHit[i],sizeof(csi_hit),1,inp);
    }
    fread(&footerVal,sizeof(uint8_t),1,inp);
    if(footerVal != 227U){
      printf("ERROR: invalid footer value in event %lu (%u)!\n", jentry, footerVal);
      exit(-1);
    }

    if(sortedEvt.header.numCsIHits>MAXNUMTIPHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIP hits (" << sortedEvt.header.numCsIHits << ")!" << endl;
      continue;
    }
    if(sortedEvt.header.numTigHits>MAXNUMTIGHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIGRESS hits (" << sortedEvt.header.numTigHits << ")!" << endl;
      continue;
    }

    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numTigHits; tipHitInd++){

      if(sortedEvt.csiHit[tipHitInd].energy > 0.0){
        int eTIP = (int)(sortedEvt.csiHit[tipHitInd].energy/MeVPerBin);
        if(eTIP>=0 && eTIP<S32K){
          mcaOut[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)][eTIP]++;
          mcaOut[NTIPRING][eTIP]++;
        }
      }
      
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  
  fclose(inp);
  
}

int main(int argc, char **argv){

  ETIP_mca_SMOL *mysort = new ETIP_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double MeVPerBin = 1.0;
  printf("Starting ETIP_mca_SMOL\n");

  if((argc != 3)&&(argc != 4)){
    cout << "Generates CsI ball mca spectra." << endl;
    cout << "Arguments: ETIP_mca_SMOL smol_file output_file MeV_per_bin" << endl;
    cout << "  *MeV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    if(argc > 3){
      MeVPerBin = atof(argv[3]);
    }
  }

  if(MeVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << MeVPerBin << ")!" << endl;
    return 0;
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << MeVPerBin << " MeV per bin." << endl;

  mysort->SortData(sfile, MeVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
