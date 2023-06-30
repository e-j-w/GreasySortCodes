//use the Makefile!

#define ComptonAngle_SEE_cxx
#include "common.h"
#include "ComptonAsym_mca_SMOL_EEGamma.h"

using namespace std;

Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position
Int_t maxECore[NTIGPOS]; //core containing the maximum energy hit, for each position
Int_t secondECore[NTIGPOS]; //core containing the 2nd highest energy hit, for each position
float maxE[NTIGPOS]; //maximum energgy deposit in keV, for each position
float mcaOut[3][S32K];

void ComptonAngle_SEE::WriteData(const char* outName){

  cout << "Writing histogram to: " << outName << endl;

  FILE *out;
  if((out = fopen(outName, "w")) == NULL){ //open the file
    cout << "ERROR: Cannot open the output file: " << outName << endl;
    return;
  }else{
    fwrite(&mcaOut,sizeof(mcaOut),1,out);
    fclose(out);
  }

}

void ComptonAngle_SEE::SortData(const char *sfile, const double ELow, const double EHigh){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  TVector3 vecBeam{0.,0.,1.}; //beam direction

  uint32_t ctsParallel = 0;
  uint32_t ctsPerp = 0;

  memset(&mcaOut,0,sizeof(mcaOut));

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    //check coincidence gate
    uint8_t pass = 0;
    uint8_t coincPos = 0;
    for(int abHitInd = 0; abHitInd < sortedEvt.header.numTigHits; abHitInd++){
      if((sortedEvt.tigHit[abHitInd].energy >= ELow)&&(sortedEvt.tigHit[abHitInd].energy <= EHigh)){
        pass = 1;
        coincPos = sortedEvt.tigHit[abHitInd].core/4;
        break;
      }
    }
    if(pass == 0){
      //no coincidence gate condition
      continue;
    }

    memset(&eABPos,0,sizeof(eABPos));
    memset(&coreFoldPos,0,sizeof(coreFoldPos));
    for(int i=0; i<NTIGPOS; i++){
      maxECore[i]=-1; //0 is a valid core number
      secondECore[i]=-1; //0 is a valid core number
      maxE[i]=0.;
    }

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].core/4 >= NTIGPOS){
        cout << "ERROR: entry " << jentry << ", hit has invalid core: " << sortedEvt.noABHit[tigHitInd].core << endl;
        exit(-1);
      }
      coreFoldPos[sortedEvt.noABHit[tigHitInd].core/4]++;
      if(sortedEvt.noABHit[tigHitInd].energy > maxE[sortedEvt.noABHit[tigHitInd].core/4]){
        if(maxECore[sortedEvt.noABHit[tigHitInd].core/4] >= 0){
          secondECore[sortedEvt.noABHit[tigHitInd].core/4]=maxECore[sortedEvt.noABHit[tigHitInd].core/4];
        }
        maxECore[sortedEvt.noABHit[tigHitInd].core/4]=sortedEvt.noABHit[tigHitInd].core;
        maxE[sortedEvt.noABHit[tigHitInd].core/4]=sortedEvt.noABHit[tigHitInd].energy;
      }else if((sortedEvt.noABHit[tigHitInd].energy > 0)&&(secondECore[sortedEvt.noABHit[tigHitInd].core/4] < 0)){
        //case where the second highest energy comes later than the highest energy
        secondECore[sortedEvt.noABHit[tigHitInd].core/4]=sortedEvt.noABHit[tigHitInd].core;
      }
      eABPos[sortedEvt.noABHit[tigHitInd].core/4] += sortedEvt.noABHit[tigHitInd].energy;
    }

    uint32_t hitPattern = 0;
    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){

      Int_t tigPos = sortedEvt.noABHit[tigHitInd].core/4;
      if(tigPos != coincPos){
        if(!(hitPattern&(1 << tigPos))){
          float eAB = eABPos[tigPos];

          if(coreFoldPos[tigPos] >= 2){
            //if((tigPos > 3)&&(tigPos < 12)){ //90 deg only
              if((eAB >= 0)&&(eAB < S32K)){
                if(maxECore[tigPos] == sortedEvt.noABHit[tigHitInd].core){
                  if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
                    TVector3 vecG1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
                    TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
                    
                    for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
                      if(secondECore[tigPos] == sortedEvt.noABHit[tigHitInd2].core){
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
                                  mcaOut[0][(int)eAB]++;
                                  mcaOut[1][(int)eAB]++;
                                }else if((angle >= -1 && angle < 20)||(angle > 160 && angle <= 181)){
                                  //cout << angle << endl;
                                  mcaOut[0][(int)eAB]--;
                                  mcaOut[2][(int)eAB]++;
                                }
                                //if((angle > 20)&&(angle < 60)){
                                  /*printf("\nPosition %i, tDiff %f\n",tigPos,tDiff);
                                  printf("Hit 1: core %2u, seg %u, vec: %.2f %.2f %.2f\n",sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg,vecG1.X(),vecG1.Y(),vecG1.Z());
                                  printf("Hit 2: core %2u, seg %u, vec: %.2f %.2f %.2f\n",sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg,vecG2.X(),vecG2.Y(),vecG2.Z());
                                  printf("Angle: %f\n",angle);*/
                                //}
                              }
                            }
                          }
                        }
                      }
                      
                    }
                  }
                }
              }
            //}
          }
          hitPattern |= (1U << tigPos); //flag this position as having been processed already
        }
      }
      
      
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;

  fclose(inp);
}

int main(int argc, char **argv){

  ComptonAngle_SEE *mysort = new ComptonAngle_SEE();

  const char *sfile, *outfile;
  printf("Starting ComptonAsym_mca_SMOL_EEGamma\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts Compton polarization asymmetry as a function of energy." << endl;
    cout << "Arguments: ComptonAsym_mca_SMOL_EEGamma smol_file EGammaLow EGammaHigh output_fmca_file" << endl;
    return 0;
  }else if(argc == 5){
    
    sfile = argv[1];
    double ELow = atof(argv[2]);
    double EHigh = atof(argv[3]);
    outfile = argv[4];

    printf("SMOL file: %s\n", sfile);
    printf("Energy gate: [%f %f] keV\n",ELow,EHigh);
    printf("Output FMCA file: %s\n", outfile);
    mysort->SortData(sfile,ELow,EHigh);
    mysort->WriteData(outfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAsym_mca_SMOL_EEGamma smol_file output_fmca_file\n");
    return 0;
  }

  return 0;
}
