//use the Makefile!

#define ComptonAngle_S_cxx
#include "common.h"
#include "ComptonAngle_SMOL_EGammaEDopp.h"

using namespace std;

#define NUM_HIST_BINS 180

int eGamma[2],eGate[2];
Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position
Double_t eABPosUnCorr[NTIGPOS]; //addback energy per position
Int_t maxECore[NTIGPOS], maxESeg[NTIGPOS];

void ComptonAngle_S::SortData(const char *sfile)
{

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  TVector3 vecBeam{0.,0.,1.}; //beam direction

  uint32_t ctsParallel = 0;
  uint32_t ctsPerp = 0;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    memset(&eABPos,0,sizeof(eABPos));
    memset(&eABPosUnCorr,0,sizeof(eABPosUnCorr));
    memset(&coreFoldPos,0,sizeof(coreFoldPos));
    memset(&maxECore,0,sizeof(maxECore));
    memset(&maxESeg,0,sizeof(maxESeg));

    Double_t maxHitE = 0.;
    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].core/4 >= NTIGPOS){
        cout << "ERROR: entry " << jentry << ", hit has invalid core: " << sortedEvt.noABHit[tigHitInd].core << endl;
        exit(-1);
      }
      coreFoldPos[sortedEvt.noABHit[tigHitInd].core/4]++;
      if(sortedEvt.noABHit[tigHitInd].energy > maxHitE){
        maxECore[sortedEvt.noABHit[tigHitInd].core/4]=sortedEvt.noABHit[tigHitInd].core;
        maxESeg[sortedEvt.noABHit[tigHitInd].core/4]=sortedEvt.noABHit[tigHitInd].seg;
        maxHitE = sortedEvt.noABHit[tigHitInd].energy;
      }
      eABPos[sortedEvt.noABHit[tigHitInd].core/4] += sortedEvt.noABHit[tigHitInd].energy;
      eABPosUnCorr[sortedEvt.noABHit[tigHitInd].core/4] += sortedEvt.noABHit[tigHitInd].energy;
    }
    for(int i=0; i<NTIGPOS; i++){
      if(eABPos[i] > MIN_TIG_EAB){
        eABPos[i] = getEDoppFusEvapManual(eABPos[i],maxECore[i],maxESeg[i],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
      }
    }

    bool passGGate = 0;
    for(int i = 0; i < NTIGPOS; i++){
      if((eABPosUnCorr[i] >= eGamma[0])&&(eABPosUnCorr[i] <= eGamma[1])){
      //if((eABPosUnCorr[i] >= 1934)&&(eABPosUnCorr[i] <= 1948)){
        passGGate = 1;
        break;
      }
    }

    if(passGGate == 1){
      for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){

        Int_t tigPos = sortedEvt.noABHit[tigHitInd].core/4;
        float eAB = eABPos[tigPos];

        if(coreFoldPos[tigPos] == 2){
          //if(!((tigPos > 3)&&(tigPos < 12))){ //non-90 deg only
          if((tigPos > 3)&&(tigPos < 12)){ //90 deg only
            if((eAB >= eGate[0])&&(eAB <= eGate[1])){

              if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
                TVector3 vecG1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
                TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
                
                for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
                  Double_t tDiff = tigHitTime(&sortedEvt,tigHitInd) - tigHitTime(&sortedEvt,tigHitInd2);
                  if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
                    //make sure both hits aren't in the same location
                    if(!(sortedEvt.noABHit[tigHitInd].core == sortedEvt.noABHit[tigHitInd2].core)){
                      //make sure both hits are in the same clover
                      if((sortedEvt.noABHit[tigHitInd2].core/4)==tigPos){
                        if(sortedEvt.noABHit[tigHitInd2].energy > MIN_TIG_EAB){
                          TVector3 vecG2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,0);
                          TVector3 norm2 = vecG2.Cross(vecG1); //norm of Compton scattering plane
                          Double_t angle = norm2.Angle(norm)*180.0/PI;
                          if((angle > 75)&&(angle < 105)){
                            ctsPerp++;
                          }else if((angle >= -1 && angle < 20)||(angle > 160 && angle <= 181)){
                            //cout << angle << endl;
                            ctsParallel++;
                          }
                          //if((angle > 20)&&(angle < 60)){
                            /*printf("\nPosition %i, tDiff %f\n",tigPos,tDiff);
                            printf("Hit 1: core %2u, seg %u, vec: %.2f %.2f %.2f\n",sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg,vecG1.X(),vecG1.Y(),vecG1.Z());
                            printf("Hit 2: core %2u, seg %u, vec: %.2f %.2f %.2f\n",sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg,vecG2.X(),vecG2.Y(),vecG2.Z());
                            printf("Angle: %f\n",angle);*/
                          //}
                          break;
                        }
                      }
                    }
                  }
                }
              }
              
              break; //only consider 1 addback hit per event
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

  cout << "Perpendicular: " << ctsPerp << endl;
  cout << "Parallel:      " << ctsParallel << endl;

  fclose(inp);
}
int main(int argc, char **argv)
{

  ComptonAngle_S *mysort = new ComptonAngle_S();

  const char *sfile;
  printf("Starting ComptonAngle_SMOL_EGammaEDopp\n");
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
    cout << "Code sorts Compton polarization asymmetry." << endl;
    cout << "Input SMOL file should contain non-addback TIGRESS data." << endl;
    cout << "Arguments: ComptonAngle_SMOL_EGammaEDopp smol_file eGammaLow eGammaHigh eGateLow eGateHigh" << endl;
    return 0;
  }else if(argc == 6){
    
    sfile = argv[1];
    eGamma[0] = atoi(argv[2]);
    eGamma[1] = atoi(argv[3]);
    eGate[0] = atoi(argv[4]);
    eGate[1] = atoi(argv[5]);

    if(eGate[0] > eGate[1]){
      //swap values
      int swapVal = eGate[1];
      eGate[1] = eGate[0];
      eGate[0] = swapVal;
    }

    printf("SMOL file: %s\n", sfile);
    cout << "Gamma coinc gate:  [" << eGamma[0] << " " << eGamma[1] << "] keV, not Doppler corrected." << endl;
    cout << "Measuring Compton angle in gate:  [" << eGate[0] << " " << eGate[1] << "] keV, Doppler corrected." << endl;
    gates = new PIDGates;
    mysort->SortData(sfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAngle_SMOL_EGammaEDopp smol_file eGammaLow eGammaHigh eGateLow eGateHigh\n");
    return 0;
  }

  return 0;
}
