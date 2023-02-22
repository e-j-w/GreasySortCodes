//use the Makefile!

#define ComptonAngle_S_cxx
#include "common.h"
#include "ComptonAngle_SMOL_EDopp.h"

using namespace std;

#define NUM_EVT_MIX 100 //number of events to scan when event mixing
#define NUM_HIST_BINS 180

int eGate[2];
Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position

void ComptonAngle_S::SortData(const char *sfile, const char *outfile)
{

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  TVector3 vecBeam{0.,0.,1.}; //beam direction

  //setup ROOT histograms
  TList *histList = new TList;
  TH1F *compAngHist = new TH1F("Compton Scattering Angle","Compton Scattering Angle",NUM_HIST_BINS,0,180);
  compAngHist->GetXaxis()->SetTitle("#theta");
  compAngHist->GetYaxis()->SetTitle("W(#theta)");
  histList->Add(compAngHist);

  uint32_t ctsParallel = 0;
  uint32_t ctsPerp = 0;

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

    memset(&eABPos,0,sizeof(eABPos));
    memset(&coreFoldPos,0,sizeof(coreFoldPos));

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){
      if(sortedEvt.tigHit[tigHitInd].core/4 >= NTIGPOS){
        cout << "ERROR: entry " << jentry << ", hit has invalid core: " << sortedEvt.tigHit[tigHitInd].core << endl;
        exit(-1);
      }
      coreFoldPos[sortedEvt.tigHit[tigHitInd].core/4]++;
      eABPos[sortedEvt.tigHit[tigHitInd].core/4] += getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
    }

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){

      Int_t tigPos = sortedEvt.tigHit[tigHitInd].core/4;
      float eAB = eABPos[tigPos];

      if(coreFoldPos[tigPos] == 2){
        //if(!((tigPos > 3)&&(tigPos < 13))){ //non-90 deg only
        //if((tigPos > 3)&&(tigPos < 13)){ //90 deg only
          if((eAB >= eGate[0])&&(eAB <= eGate[1])){

            Double_t eDopp = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);

            if(eDopp > 30.){
              TVector3 vecG1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
              TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
              
              for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
                Double_t tDiff = tigHitTime(&sortedEvt,tigHitInd) - tigHitTime(&sortedEvt,tigHitInd2);
                if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
                  //make sure both hits aren't in the same location
                  if(!(sortedEvt.tigHit[tigHitInd].core == sortedEvt.tigHit[tigHitInd2].core)){
                    //make sure both hits are in the same clover
                    if((sortedEvt.tigHit[tigHitInd2].core/4)==tigPos){
                      Double_t eDopp2 = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
                      if(eDopp2 > 30.){
                        TVector3 vecG2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
                        TVector3 norm2 = vecG2.Cross(vecG1); //norm of Compton scattering plane
                        Double_t angle = norm2.Angle(norm)*180/PI;
                        compAngHist->Fill(angle);
                        if((angle > 75)&&(angle < 105)){
                          ctsPerp++;
                        }else if((angle >= 0 && angle < 20)||(angle > 160 && angle <= 180)){
                          //cout << angle << endl;
                          ctsParallel++;
                        }
                        /*if((angle > 20)&&(angle < 60)){
                          printf("\nPosition %i, tDiff %f\n",tigPos,tDiff);
                          printf("Hit 1: core %u, seg %u, energy %f\n",sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg,eDopp);
                          printf("Hit 2: core %u, seg %u, energy %f\n",sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg,eDopp2);
                          //printf("Vec: [%f %f %f]\n",vecCS.X(),vecCS.Y(),vecCS.Z());
                          printf("Angle: %f\n",angle);
                        }*/
                        break;
                      }
                    }
                  }
                }
              }
            }
            
            break; //only consider 1 addback hit per event
          }
        //}
      }
      
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;

  //dump hists to arrays
  cout << "Angle   Cts" << endl;
  for(int i=0; i<=compAngHist->GetNbinsX(); i++){
    if(compAngHist->GetBinContent(i) > 0.){
      cout << compAngHist->GetBinCenter(i) << " " << compAngHist->GetBinContent(i) << endl;
    }
  }

  cout << "Perpendicular: " << ctsPerp << endl;
  cout << "Parallel:      " << ctsParallel << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  histList->Write();
  myfile->Write();
  myfile->Close();

  fclose(inp);
}
int main(int argc, char **argv)
{

  ComptonAngle_S *mysort = new ComptonAngle_S();

  const char *sfile, *outfile;
  printf("Starting ComptonAngle_SMOL_EDopp\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts DCO ratios." << endl;
    cout << "Arguments: ComptonAngle_SMOL_EDopp smol_file eGateLow eGateHigh output_root_file" << endl;
    return 0;
  }else if(argc == 5){
    
    sfile = argv[1];
    eGate[0] = atoi(argv[2]);
    eGate[1] = atoi(argv[3]);
    outfile = argv[4];

    if(eGate[0] > eGate[1]){
      //swap values
      int swapVal = eGate[1];
      eGate[1] = eGate[0];
      eGate[0] = swapVal;
    }

    printf("SMOL file: %s\n", sfile);
    cout << "Gate:  [" << eGate[0] << " " << eGate[1] << "] keV" << endl;
    printf("Output file: %s\n", outfile);
    gates = new PIDGates;
    mysort->SortData(sfile,outfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAngle_SMOL_EDopp smol_file eGateLow eGateHigh output_root_file\n");
    return 0;
  }

  return 0;
}
