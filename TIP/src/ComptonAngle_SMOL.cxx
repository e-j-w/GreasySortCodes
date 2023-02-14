//use the Makefile!

#define ComptonAngle_S_cxx
#include "common.h"
#include "ComptonAngle_SMOL.h"

using namespace std;

#define NUM_EVT_MIX 100 //number of events to scan when event mixing
#define NUM_HIST_BINS 180

int eGate[2];
Int_t coreFoldPos[NTIGPOS]; //number of core hits per position
Double_t eABPos[NTIGPOS]; //addback energy per position
sorted_evt pastEvt[NUM_EVT_MIX];

double xRaw[NUM_HIST_BINS], yRaw[NUM_HIST_BINS];
double xEvtMix[NUM_HIST_BINS], yEvtMix[NUM_HIST_BINS];
double xNorm[NUM_HIST_BINS], yNorm[NUM_HIST_BINS], eNorm[NUM_HIST_BINS];

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
  TH1F *compAngHistEvtMix = new TH1F("Event Mixed Compton Scattering Angle","Event Mixed Compton Scattering Angle",NUM_HIST_BINS,0,180);
  compAngHistEvtMix->GetXaxis()->SetTitle("#theta");
  compAngHistEvtMix->GetYaxis()->SetTitle("W(#theta)");
  histList->Add(compAngHistEvtMix);

  printf("\nSorting events...\n");
  Int_t evtMixCtr = 0;
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    memset(&sortedEvt,0,sizeof(sorted_evt));
    footerVal = 0;
    fread(&sortedEvt.header,sizeof(evt_header),1,inp);
    for(int i = 0; i<sortedEvt.header.numTigABHits;i++){
      fread(&sortedEvt.tigHit[i],sizeof(tigab_hit),1,inp);
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
    if(sortedEvt.header.numTigABHits>MAXNUMTIGHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIGRESS hits (" << sortedEvt.header.numTigABHits << ")!" << endl;
      continue;
    }

    memcpy(&pastEvt[evtMixCtr],&sortedEvt,sizeof(sorted_evt)); //copy event to past event array
    memset(&eABPos,0,sizeof(eABPos));
    memset(&coreFoldPos,0,sizeof(coreFoldPos));

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigABHits; tigHitInd++){
      if(sortedEvt.tigHit[tigHitInd].core/4 >= NTIGPOS){
        cout << "ERROR: entry " << jentry << ", hit has invalid core: " << sortedEvt.tigHit[tigHitInd].core << endl;
        exit(-1);
      }
      coreFoldPos[sortedEvt.tigHit[tigHitInd].core/4]++;
      eABPos[sortedEvt.tigHit[tigHitInd].core/4] += sortedEvt.tigHit[tigHitInd].energy;
    }

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigABHits; tigHitInd++){

      Int_t tigPos = sortedEvt.tigHit[tigHitInd].core/4;
      float eAB = eABPos[tigPos];

      if(coreFoldPos[tigPos] == 2){
        if((tigPos > 3)&&(tigPos < 13)){
          if(eAB > MIN_TIG_EAB){
            if((eAB >= eGate[0])&&(eAB <= eGate[1])){

              if(sortedEvt.tigHit[tigHitInd].energy > MIN_TIG_EAB){
                TVector3 vecG1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
                TVector3 norm = vecG1.Cross(vecBeam); //norm of reaction plane
                
                for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigABHits; tigHitInd2++){
                  Double_t tDiff = sortedEvt.tigHit[tigHitInd].timeNs - sortedEvt.tigHit[tigHitInd2].timeNs;
                  if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
                    //make sure both hits aren't in the same location
                    if(!(sortedEvt.tigHit[tigHitInd].core == sortedEvt.tigHit[tigHitInd2].core)){
                      //make sure both hits are in the same clover
                      if((sortedEvt.tigHit[tigHitInd2].core/4)==tigPos){
                        if(sortedEvt.tigHit[tigHitInd2].energy > MIN_TIG_EAB){
                          TVector3 vecG2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
                          TVector3 norm2 = vecG2.Cross(vecG1); //norm of Compton scattering plane
                          Double_t angle = norm2.Angle(norm)*180/PI;
                          compAngHist->Fill(angle);
                          /*if((angle > 20)&&(angle < 60)){
                            printf("\nPosition %i, tDiff %f\n",tigPos,tDiff);
                            printf("Hit 1: core %u, seg %u, energy %f\n",sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg,sortedEvt.tigHit[tigHitInd].energy);
                            printf("Hit 2: core %u, seg %u, energy %f\n",sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg,sortedEvt.tigHit[tigHitInd2].energy);
                            //printf("Vec: [%f %f %f]\n",vecCS.X(),vecCS.Y(),vecCS.Z());
                            printf("Angle: %f\n",angle);
                          }*/
                          break;
                        }
                      }
                    }
                  }
                }

                //event mixing
                for(int evtMixInd = 0; evtMixInd < NUM_EVT_MIX; evtMixInd++){
                  //handle early events where event mix array isn't full yet
                  if(evtMixInd < jentry){
                    //don't mix event with itself
                    if(evtMixInd != evtMixCtr){
                      for(int tigHitInd2 = 0; tigHitInd2 < pastEvt[evtMixInd].header.numTigABHits; tigHitInd2++){
                        //make sure both hits aren't in the same location
                        if(!(sortedEvt.tigHit[tigHitInd].core == pastEvt[evtMixInd].tigHit[tigHitInd2].core)){
                          //make sure both hits are in the same clover
                          if((pastEvt[evtMixInd].tigHit[tigHitInd2].core/4)==tigPos){
                            if(pastEvt[evtMixInd].tigHit[tigHitInd2].energy > MIN_TIG_EAB){
                              TVector3 vecG2 = getTigVector(pastEvt[evtMixInd].tigHit[tigHitInd2].core,0);
                              TVector3 norm2 = vecG2.Cross(vecG1); //norm of Compton scattering plane
                              Double_t angle = norm2.Angle(norm);
                              compAngHistEvtMix->Fill(angle*180/PI);
                            }
                          }
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

    //increment event mix array index
    evtMixCtr++;
    if(evtMixCtr >= NUM_EVT_MIX){
      evtMixCtr = 0;
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;

  //dump hists to arrays
  int bin = 0;
  for(int i=0; i<compAngHist->GetNbinsX(); i++){
    if(compAngHist->GetBinContent(i) > 0.){
      if(bin<NUM_HIST_BINS){
        xRaw[bin]=compAngHist->GetBinCenter(i);
        yRaw[bin]=compAngHist->GetBinContent(i);
        bin++;
      }
    }
    
  }
  bin = 0;
  for(int i=0; i<compAngHistEvtMix->GetNbinsX(); i++){
    if(compAngHist->GetBinContent(i) > 0.){
      if(bin<NUM_HIST_BINS){
        xEvtMix[bin]=compAngHistEvtMix->GetBinCenter(i);
        yEvtMix[bin]=compAngHistEvtMix->GetBinContent(i);
      }
      bin++;
    }
  }
  //generate normalized data
  for(int i=0; i<bin; i++){
    if(i<NUM_HIST_BINS){
      xNorm[i]=xRaw[i];
      if((yEvtMix[i] > 0.)&&(yRaw[i] > 0.)){
        yNorm[i] = (yRaw[i]/yEvtMix[i]);
        eNorm[i] = yNorm[i]*sqrt((1./yRaw[i]) + (1./yEvtMix[i]));
        cout << xNorm[i] << " " << yRaw[i] << " " << yEvtMix[i] << " " << yNorm[i] << endl;
      }else{
        yNorm[i] = 0.;
        eNorm[i] = 0.;
      }
    }
  }

  //plot data with error bars 
	TGraphErrors *grNorm = new TGraphErrors(bin,xNorm,yNorm,nullptr,eNorm);
	grNorm->SetTitle("Normalized Compton scattering angular distribution");
	grNorm->GetYaxis()->SetTitle("W(#theta)");
	grNorm->GetYaxis()->SetTitleSize(0.085);
	grNorm->GetYaxis()->SetTitleOffset(0.5);
	grNorm->GetYaxis()->CenterTitle();
	grNorm->GetYaxis()->SetLabelSize(0.065);
	grNorm->GetYaxis()->SetNdivisions(6);
	grNorm->GetYaxis()->SetDecimals();
	//grNorm->GetYaxis()->SetRangeUser(0.4,1.6);
	grNorm->GetXaxis()->SetTitle("#theta");
	grNorm->GetXaxis()->SetTitleSize(0.085);
	grNorm->GetXaxis()->SetTitleOffset(0.85);
	grNorm->GetXaxis()->CenterTitle();
	grNorm->GetXaxis()->SetLabelSize(0.065);
	grNorm->GetXaxis()->SetNdivisions(5);
	grNorm->GetXaxis()->SetLimits(0,180);
	grNorm->SetMarkerColor(1);
	grNorm->SetMarkerStyle(21);
	grNorm->SetMarkerSize(1.5);
  histList->Add(grNorm);

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
  printf("Starting ComptonAngle_SMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts DCO ratios." << endl;
    cout << "Arguments: ComptonAngle_SMOL smol_file eGateLow eGateHigh output_root_file" << endl;
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
    mysort->SortData(sfile,outfile);

  }else{
    printf("ERROR: Improper number of arguments!\nArguments: ComptonAngle_SMOL smol_file eGateLow eGateHigh output_root_file\n");
    return 0;
  }

  return 0;
}
