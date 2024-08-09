//use the Makefile!

#define RDCO_S_cxx
#include "common.h"
#include "RDCO_SMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];
int eGate[2], eGamma[2];

Int_t ctsRingGateRing[NTIGSEGRING-1][NTIGSEGRING-1];

void RDCO_S::SortData(char const *sfile)
{

  memset(ctsRingGateRing,0,sizeof(ctsRingGateRing));

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

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){
      float addE1 = sortedEvt.tigHit[tigHitInd].energy;

      if(addE1 > MIN_TIG_EAB){
        if((addE1 >= eGate[0])&&(addE1 <= eGate[1])){
          TVector3 vecGate = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
          Int_t ringGate = getTIGRESSSegmentRing(vecGate.Theta()*180.0/PI);
          for(int tigHitInd2 = 0; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            if(tigHitInd2 != tigHitInd){
              float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
              if(addE2 > MIN_TIG_EAB){
                if((addE2 >= eGamma[0])&&(addE2 <= eGamma[1])){
                  //cascade identified
                  TVector3 vecGamma = getTigVector(sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg);
                  Int_t ringGamma = getTIGRESSSegmentRing(vecGamma.Theta()*180.0/PI);
                  if((ringGate < (NTIGSEGRING-1))&&(ringGamma < (NTIGSEGRING-1))){
                    ctsRingGateRing[ringGamma][ringGate]++;
                  }
                }
              }
            }
          }
        }
      }
      
    }
    

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;

  Double_t RDCONum = 0.;
  Double_t RDCODen = 0.;

  for(Int_t i=0; i<NTIGSEGRING-1; i++){
    for(Int_t j=0; j<NTIGSEGRING-1; j++){
      if(i!=j){
        if(i<4 || i>7){ //exclude 90 deg
          if(j>=4 && j<=7){ //only 90 deg
            RDCONum += ctsRingGateRing[i][j];
            RDCODen += ctsRingGateRing[j][i];
          }
        }
      }
    }
  }

  Double_t RDCO = RDCONum/RDCODen;
  cout << "Counts: " << RDCONum << " " << RDCODen << endl;
  cout << "RDCO: " << RDCO << " +/- " << RDCO*sqrt((1./RDCONum) + (1./RDCODen)) << endl;

  fclose(inp);
}
int main(int argc, char **argv)
{

  RDCO_S *mysort = new RDCO_S();

  const char *sfile;
  printf("Starting RDCO_SMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts DCO ratios." << endl;
    cout << "Arguments: RDCO_SMOL smol_file eGateLow eGateHigh eGammaLow eGammaHigh" << endl;
    return 0;
  }else if(argc == 6){
    sfile = argv[1];
    eGate[0] = atoi(argv[2]);
    eGate[1] = atoi(argv[3]);
    eGamma[0] = atoi(argv[4]);
    eGamma[1] = atoi(argv[5]);
    printf("SMOL file: %s\n", sfile); 
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: RDCO_SMOL smol_file eGateLow eGateHigh eGammaLow eGammaHigh\n");
    return 0;
  }

  if(eGate[0] > eGate[1]){
    //swap values
    int swapVal = eGate[1];
    eGate[1] = eGate[0];
    eGate[0] = swapVal;
  }
  if(eGamma[0] > eGamma[1]){
    //swap values
    int swapVal = eGamma[1];
    eGamma[1] = eGamma[0];
    eGamma[0] = swapVal;
  }

  cout << "Gate:  [" << eGate[0] << " " << eGate[1] << "] keV" << endl;
  cout << "Gamma: [" << eGamma[0] << " " << eGamma[1] << "] keV" << endl;

  mysort->SortData(sfile);

  return 0;
}
