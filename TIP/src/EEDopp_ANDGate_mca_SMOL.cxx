//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define EEDopp_ANDGate_mca_SMOL_cxx
#include "common.h"
#include "EEDopp_ANDGate_mca_SMOL.h"

using namespace std;

void EEDopp_ANDGate_mca_SMOL::WriteData(const char* outName){

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

void EEDopp_ANDGate_mca_SMOL::SortData(char const *sfile, const uint8_t numEGates, const double eLow[10], const double eHigh[10], const double keVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    uint8_t hitsInGates = 0; //hitpattern for flagging hits in one of the energy gates
    uint8_t gatesHit = 0; //flags of which energy gates were satisfied

    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){
      int eDopp = (int)(getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates));
      for(uint8_t i=0; i<numEGates;i++){
        if(!(gatesHit & (1 << i))){
          if((eDopp >= eLow[i])&&(eDopp <= eHigh[i])){
            gatesHit |= (1 << i);
            hitsInGates |= (1 << tigHitIndAB);
            break;
          }
        }
      }
    }

    uint8_t numGatesHit = 0;
    for(uint8_t i=0; i<numEGates;i++){
      if(gatesHit & (1 << i)){
        numGatesHit++;
      }
    }
    if(numGatesHit == numEGates){
      for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){
        if(!(hitsInGates & (1 << tigHitIndAB))){
          //hit not in any gates, add to spectrum
          int eDopp = (int)(getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates));
          if(eDopp>=0 && eDopp<S32K){
            double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180./PI;
            mcaOut[getTIGRESSRing(theta)+1][eDopp]++;
            mcaOut[getTIGRESSSegmentRing(theta)+7][eDopp]++;
            mcaOut[0][eDopp]++;
          }
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

  EEDopp_ANDGate_mca_SMOL *mysort = new EEDopp_ANDGate_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double keVPerBin = 0.0;
  double eLow[10], eHigh[10];
  uint8_t numEGates = 0;
  printf("Starting EEDopp_ANDGate_mca_SMOL\n");

  if((argc < 6)||((argc % 2)!=0)){
    cout << "Generates TIGRESS mca spectra for PID and time separated data." << endl;
    cout << "Arguments: EEDopp_ANDGate_mca_SMOL smol_file output_file keV_per_bin EGateLow1 EGateHigh1 EGateLow2 EGateHigh2 ..." << endl;
    return 0;
  }else if(argc > 24){
    cout << "ERROR: Too many parameters, maximum 10 energy gates allowed." << endl;
    return 0;
  }else{
    numEGates = (argc - 4)/2;
    sfile = argv[1];
    outfile = argv[2];
    keVPerBin = atof(argv[3]);
    if(keVPerBin <= 0.0){
      cout << "ERROR: Invalid keV/bin factor (" << keVPerBin << ")!" << endl;
      return 0;
    }
    for(uint8_t i=0;i<numEGates;i++){
      eLow[i] = atof(argv[4+2*i]);
      eHigh[i] = atof(argv[5+2*i]);
      if(eLow[i] > eHigh[i]){
        //swap energy gate bounds
        double tmp = eHigh[i];
        eHigh[i] = eLow[i];
        eLow[i] = tmp;
      }
      cout << "Energy gate " << (int)i << ": [" << eLow[i] << " " << eHigh[i] << "] keV." << endl;
    }
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum
  gates = new PIDGates;

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << keVPerBin << " keV per bin." << endl;

  mysort->SortData(sfile, numEGates, eLow, eHigh, keVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
