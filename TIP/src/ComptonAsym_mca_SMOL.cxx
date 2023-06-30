//use the Makefile!

#define ComptonAngle_S_cxx
#include "common.h"
#include "ComptonAsym_mca_SMOL.h"

using namespace std;

Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position
Int_t maxECore[NTIGPOS]; //core containing the maximum energy hit, for each position
Int_t secondECore[NTIGPOS]; //core containing the 2nd highest energy hit, for each position
float maxE[NTIGPOS]; //maximum energgy deposit in keV, for each position
float mcaOut[3][S32K];

void ComptonAngle_S::WriteData(const char* outName){

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

  memset(&mcaOut,0,sizeof(mcaOut));

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
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
      if(!(hitPattern&(1 << tigPos))){
        float eAB = eABPos[tigPos];

        if(coreFoldPos[tigPos] >= 2){
          //if((tigPos < 4)||(tigPos > 11)){
          //if((tigPos > 3)&&(tigPos < 12)){ //90 deg only
            if((eAB >= 0)&&(eAB < S32K)){
              if(maxECore[tigPos] == sortedEvt.noABHit[tigHitInd].core){
                if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
                  TVector3 vecG1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg);
                  TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
                  
                  for(int tigHitInd2 = 0; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
                    if(secondECore[tigPos] == sortedEvt.noABHit[tigHitInd2].core){

                      /*//check segments
                      uint8_t pass = 0;
                      if((sortedEvt.noABHit[tigHitInd].seg < 5)&&(sortedEvt.noABHit[tigHitInd2].seg < 5)){
                        pass = 1;
                      }else if((sortedEvt.noABHit[tigHitInd].seg > 4)&&(sortedEvt.noABHit[tigHitInd2].seg > 4)){
                        pass = 1;
                      }
                      if(pass == 0){
                        continue;
                      }*/

                      Double_t tDiff = noABHitTime(&sortedEvt,tigHitInd) - noABHitTime(&sortedEvt,tigHitInd2);
                      if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
                        //make sure both hits aren't in the same location
                        if(!(sortedEvt.noABHit[tigHitInd].core == sortedEvt.noABHit[tigHitInd2].core)){
                          //make sure both hits are in the same clover
                          if((sortedEvt.noABHit[tigHitInd2].core/4)==tigPos){
                            if(sortedEvt.noABHit[tigHitInd2].energy > MIN_TIG_EAB){
                              TVector3 vecG2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg);
                              TVector3 norm2 = vecG2.Cross(vecG1); //norm of Compton scattering plane
                              Double_t angle = norm2.Angle(norm)*180.0/PI;
                              if((angle > 85)&&(angle < 95)){
                                mcaOut[0][(int)eAB]++;
                                mcaOut[1][(int)eAB]++;
                              }else if((angle >= -1 && angle < 5)||(angle > 175 && angle <= 181)){
                                //cout << angle << endl;
                                mcaOut[0][(int)eAB]--;
                                mcaOut[2][(int)eAB]++;
                              }
                              //if((angle > 20)&&(angle < 60)){
                              /*if(eAB>3555 && eAB<3570){
                                printf("\nPosition %i, tDiff %f, addback E %f\n",tigPos,tDiff,eAB);
                                printf("Hit 1 (index %i): energy %8.3f, core %2u, seg %u, vec: %.2f %.2f %.2f\n",tigHitInd,sortedEvt.noABHit[tigHitInd].energy,sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg,vecG1.X(),vecG1.Y(),vecG1.Z());
                                printf("Hit 2 (index %i): energy %8.3f, core %2u, seg %u, vec: %.2f %.2f %.2f\n",tigHitInd2,sortedEvt.noABHit[tigHitInd2].energy,sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg,vecG2.X(),vecG2.Y(),vecG2.Z());
                                printf("Dist: %f\n",(vecG1 - vecG2).Mag());
                                printf("Angle: %f\n",angle);
                                if((vecG1 - vecG2).Mag() > 118){
                                  getc(stdin);
                                }

                                printf("Non-addback hits:\n");
                                for (int i=0; i<sortedEvt.header.numNoABHits; i++){
                                  printf("Hit %i: energy %8.3f, core %2u\n",i,sortedEvt.noABHit[i].energy,sortedEvt.noABHit[i].core);
                                }
                                printf("Addback hits:\n");
                                for (int i=0; i<sortedEvt.header.numNoABHits; i++){
                                  printf("Hit %i: energy %8.3f, core %2u\n",i,sortedEvt.noABHit[i].energy,sortedEvt.noABHit[i].core);
                                }
                                printf("\n");
                              }*/
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

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;

  fclose(inp);
}

int main(int argc, char **argv){

  ComptonAngle_S *mysort = new ComptonAngle_S();

  const char *sfile, *outfile;
  printf("Starting ComptonAsym_mca_SMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts Compton polarization asymmetry as a function of energy." << endl;
    cout << "Arguments: ComptonAsym_mca_SMOL smol_file output_fmca_file" << endl;
    return 0;
  }else if(argc == 3){
    
    sfile = argv[1];
    outfile = argv[2];

    printf("SMOL file: %s\n", sfile);
    printf("Output FMCA file: %s\n", outfile);
    mysort->SortData(sfile);
    mysort->WriteData(outfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAsym_mca_SMOL smol_file output_fmca_file\n");
    return 0;
  }

  return 0;
}
