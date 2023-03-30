//use the Makefile!

#define ComptonAngle_S_cxx
#include "common.h"
#include "ComptonAngle_SMOL.h"

using namespace std;

#define NUM_HIST_BINS 180

int eGate[2];
Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position

void ComptonAngle_S::SortData(const char *sfile){

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
    memset(&coreFoldPos,0,sizeof(coreFoldPos));

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].core/4 >= NTIGPOS){
        cout << "ERROR: entry " << jentry << ", hit has invalid core: " << sortedEvt.noABHit[tigHitInd].core << endl;
        exit(-1);
      }
      coreFoldPos[sortedEvt.noABHit[tigHitInd].core/4]++;
      eABPos[sortedEvt.noABHit[tigHitInd].core/4] += sortedEvt.noABHit[tigHitInd].energy;
    }

    uint32_t hitPattern = 0;
    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){

      Int_t tigPos = sortedEvt.noABHit[tigHitInd].core/4;
      if(!(hitPattern&(1 << tigPos))){
        float eAB = eABPos[tigPos];

        if(coreFoldPos[tigPos] == 2){
          //if((tigPos > 3)&&(tigPos < 12)){ //90 deg only
            if((eAB >= eGate[0])&&(eAB <= eGate[1])){

              if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
                TVector3 vecG1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
                TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
                
                for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
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
              
              hitPattern |= (1U << tigPos); //flag this position as having been processed already
            }
          //}
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
  printf("Starting ComptonAngle_SMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts Compton polarization asymmetry." << endl;
    cout << "Arguments: ComptonAngle_SMOL smol_file eGateLow eGateHigh" << endl;
    return 0;
  }else if(argc == 4){
    
    sfile = argv[1];
    eGate[0] = atoi(argv[2]);
    eGate[1] = atoi(argv[3]);

    if(eGate[0] > eGate[1]){
      //swap values
      int swapVal = eGate[1];
      eGate[1] = eGate[0];
      eGate[0] = swapVal;
    }

    printf("SMOL file: %s\n", sfile);
    cout << "Gate:  [" << eGate[0] << " " << eGate[1] << "] keV" << endl;
    mysort->SortData(sfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAngle_SMOL smol_file eGateLow eGateHigh\n");
    return 0;
  }

  return 0;
}