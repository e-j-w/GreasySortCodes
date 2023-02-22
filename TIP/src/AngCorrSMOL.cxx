//use the Makefile!

#define AngCorrS_cxx
#include "common.h"
#include "AngCorrSMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];
int eGamma1[2], eGamma2[2], eBG1[2], eBG2[2];
double xCore[NUM_HIST_BINS], yCoreRaw[NUM_HIST_BINS], yCoreBG[NUM_HIST_BINS], yCoreEvtMix[NUM_HIST_BINS];
double eCoreRaw[NUM_HIST_BINS], eCoreBG[NUM_HIST_BINS], eCoreEvtMix[NUM_HIST_BINS];
double yCoreNorm[NUM_HIST_BINS], eCoreNorm[NUM_HIST_BINS];


void AngCorr::SortData(char const *sfile, char const *outfile)
{
  Initialise();

  double bgScaling = ((eGamma1[1]-eGamma1[0])/(1.0*(eBG1[1]-eBG1[0])))*((eGamma2[1]-eGamma2[0])/(1.0*(eBG2[1]-eBG2[0])));
  //bgScaling = 0;
  cout << "Background scaling factor: " << bgScaling << endl;

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  //for event mixing
  tig_hit pastHit[NUM_PAST_HIT];
  memset(pastHit,0,sizeof(pastHit));
  Long64_t pastHitEvt[NUM_PAST_HIT];
  memset(pastHitEvt,0,sizeof(pastHitEvt));
  int pastEvtInd = 0;

  /*cout << "TIGRESS positions: " << endl;
  for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
      TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
      cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
  }*/

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

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){
      float addE1 = sortedEvt.tigHit[tigHitInd].energy;

      if(addE1 > MIN_TIG_EAB){
        if((addE1 >= eGamma1[0])&&(addE1 <= eGamma1[1])){
          for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
            if((addE2 >= eGamma2[0])&&(addE2 <= eGamma2[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrRaw->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreRaw->Fill(cos(coreangle));
            }
          }
          memcpy(&pastHit[pastEvtInd],&sortedEvt.tigHit[tigHitInd],sizeof(tig_hit));
          pastHitEvt[pastEvtInd] = jentry;
          pastEvtInd++;
          if(pastEvtInd >= NUM_PAST_HIT){
            pastEvtInd = 0; //roll over
          }
        }else if((addE1 >= eGamma2[0])&&(addE1 <= eGamma2[1])){
          for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
            if((addE2 >= eGamma1[0])&&(addE2 <= eGamma1[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrRaw->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreRaw->Fill(cos(coreangle));
            }
          }
          memcpy(&pastHit[pastEvtInd],&sortedEvt.tigHit[tigHitInd],sizeof(tig_hit));
          pastHitEvt[pastEvtInd] = jentry;
          pastEvtInd++;
          if(pastEvtInd >= NUM_PAST_HIT){
            pastEvtInd = 0; //roll over
          }
        }
        if((addE1 >= eBG1[0])&&(addE1 <= eBG1[1])){
          for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
            if((addE2 >= eBG2[0])&&(addE2 <= eBG2[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrBG->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreBG->Fill(cos(coreangle));
            }
          }
          memcpy(&pastHit[pastEvtInd],&sortedEvt.tigHit[tigHitInd],sizeof(tig_hit));
          pastHitEvt[pastEvtInd] = jentry;
          pastEvtInd++;
          if(pastEvtInd >= NUM_PAST_HIT){
            pastEvtInd = 0; //roll over
          }
        }else if((addE1 >= eBG2[0])&&(addE1 <= eBG2[1])){
          for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
            if((addE2 >= eBG1[0])&&(addE2 <= eBG1[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,sortedEvt.tigHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrBG->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.tigHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreBG->Fill(cos(coreangle));
            }
          }
          memcpy(&pastHit[pastEvtInd],&sortedEvt.tigHit[tigHitInd],sizeof(tig_hit));
          pastHitEvt[pastEvtInd] = jentry;
          pastEvtInd++;
          if(pastEvtInd >= NUM_PAST_HIT){
            pastEvtInd = 0; //roll over
          }
        }

        //handle event mixing
        if((addE1 >= eGamma1[0])&&(addE1 <= eGamma1[1])){
          for(int tigHitInd2 = 0; tigHitInd2 < NUM_PAST_HIT; tigHitInd2++){
            if((pastHitEvt[tigHitInd2] < jentry)&&(pastHitEvt[tigHitInd2] > (jentry-1000))){
              //hit is from a different event
              float addE2 = pastHit[tigHitInd2].energy;
              if((addE2 >= eGamma2[0])&&(addE2 <= eGamma2[1])){
                //event mixed cascade identified
                TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
                TVector3 vec2 = getTigVector(pastHit[tigHitInd2].core,pastHit[tigHitInd2].seg);
                Double_t angle = vec1.Angle(vec2);
                angCorrEvtMix->Fill(cos(angle));
                TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
                TVector3 corevec2 = getTigVector(pastHit[tigHitInd2].core,0);
                Double_t coreangle = corevec1.Angle(corevec2);
                angCorrCoreEvtMix->Fill(cos(coreangle));
                evtMixDist->Fill(jentry - pastHitEvt[tigHitInd2]);
              }
            }
          }
        }else if((addE1 >= eGamma2[0])&&(addE1 <= eGamma2[1])){
          for(int tigHitInd2 = 0; tigHitInd2 < NUM_PAST_HIT; tigHitInd2++){
            if((pastHitEvt[tigHitInd2] < jentry)&&(pastHitEvt[tigHitInd2] > (jentry-1000))){
              //hit is from a different event
              float addE2 = pastHit[tigHitInd2].energy;
              if((addE2 >= eGamma1[0])&&(addE2 <= eGamma1[1])){
                //event mixed cascade identified
                TVector3 vec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
                TVector3 vec2 = getTigVector(pastHit[tigHitInd2].core,pastHit[tigHitInd2].seg);
                Double_t angle = vec1.Angle(vec2);
                angCorrEvtMix->Fill(cos(angle));
                TVector3 corevec1 = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
                TVector3 corevec2 = getTigVector(pastHit[tigHitInd2].core,0);
                Double_t coreangle = corevec1.Angle(corevec2);
                angCorrCoreEvtMix->Fill(cos(coreangle));
                evtMixDist->Fill(jentry - pastHitEvt[tigHitInd2]);
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

  //dump histograms to arrays
  const int contractFactor = 2;
  int contractCounter = 0;
  int bin=0, numBins=0;
  xCore[bin]=0.;
  yCoreRaw[bin]=0.;
  eCoreRaw[bin]=0.;
  yCoreBG[bin]=0.;
  eCoreBG[bin]=0.;
  yCoreEvtMix[bin]=0.;
  eCoreEvtMix[bin]=0.;
  for(int i=0; i<angCorrCoreRaw->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrCoreRaw->GetBinContent(i) > 0.){
        //cout << angCorrCoreRaw->GetBinCenter(i) << " " << angCorrCoreRaw->GetBinContent(i) << endl;
        xCore[bin]+=angCorrCoreRaw->GetBinCenter(i);
        yCoreRaw[bin]+=angCorrCoreRaw->GetBinContent(i);
        contractCounter++;
        if(contractCounter >= contractFactor){
          xCore[bin] /= contractFactor;
          //cout << "contracted to " <<  xCore[bin] << " " << yCoreRaw[bin] << endl;
          eCoreRaw[bin] = sqrt(yCoreRaw[bin]);
          bin++;
          contractCounter = 0;
          xCore[bin]=0.;
          yCoreRaw[bin]=0.;
          eCoreRaw[bin]=0.;
        }
      }
    }
  }
  numBins=bin;
  bin=0;
  contractCounter = 0;
  for(int i=0; i<angCorrCoreBG->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrCoreRaw->GetBinContent(i) > 0.){
        if(angCorrCoreBG->GetBinCenter(i) != angCorrCoreRaw->GetBinCenter(i)){
          printf("x-value mismatch between raw and BG at bin %i: %f %f\n",i,xCore[bin],angCorrCoreBG->GetBinCenter(i));
        }
        yCoreBG[bin]+=angCorrCoreBG->GetBinContent(i);
        contractCounter++;
        if(contractCounter >= contractFactor){
          eCoreBG[bin] = sqrt(yCoreBG[bin]);
          bin++;
          contractCounter = 0;
          yCoreBG[bin]=0.;
          eCoreBG[bin]=0.;
        }
      }
    }
  }
  bin=0;
  contractCounter = 0;
  for(int i=0; i<angCorrCoreEvtMix->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrCoreEvtMix->GetBinContent(i) > 0.){
        if(angCorrCoreEvtMix->GetBinCenter(i) != angCorrCoreRaw->GetBinCenter(i)){
          printf("x-value mismatch between raw and event mixed at bin %i: %f %f\n",i,xCore[bin],angCorrCoreEvtMix->GetBinCenter(i));
        }
        yCoreEvtMix[bin]+=angCorrCoreEvtMix->GetBinContent(i);
        contractCounter++;
        if(contractCounter >= contractFactor){
          eCoreEvtMix[bin] = sqrt(yCoreEvtMix[bin]);
          bin++;
          contractCounter = 0;
          yCoreEvtMix[bin]=0.;
          eCoreEvtMix[bin]=0.;
        }
      }
    }
  }

  //generate normalized data
  double avg = 0.;
  for(int i=0; i<numBins; i++){
    yCoreNorm[i] = (yCoreRaw[i] - bgScaling*yCoreBG[i])/yCoreEvtMix[i];
    eCoreNorm[i] = ((sqrt(pow(eCoreRaw[i],2) + pow(bgScaling*eCoreBG[i],2)))/(yCoreRaw[i] - bgScaling*yCoreBG[i]));
    eCoreNorm[i] += eCoreEvtMix[i]/yCoreEvtMix[i];
    eCoreNorm[i] *= yCoreNorm[i];
    avg += yCoreNorm[i];
  }
  avg /= numBins;
  for(int i=0; i<numBins; i++){
    yCoreNorm[i] /= avg;
    eCoreNorm[i] /= avg;
  }

  //plot data with error bars 
	TGraphErrors *grCoreNorm = new TGraphErrors(numBins-3,xCore,yCoreNorm,nullptr,eCoreNorm);
	grCoreNorm->SetTitle("Normalized angular correlation");
	grCoreNorm->GetYaxis()->SetTitle("W(#theta)");
	grCoreNorm->GetYaxis()->SetTitleSize(0.085);
	grCoreNorm->GetYaxis()->SetTitleOffset(0.5);
	grCoreNorm->GetYaxis()->CenterTitle();
	grCoreNorm->GetYaxis()->SetLabelSize(0.065);
	grCoreNorm->GetYaxis()->SetNdivisions(6);
	grCoreNorm->GetYaxis()->SetDecimals();
	grCoreNorm->GetYaxis()->SetRangeUser(0.4,1.6);
	grCoreNorm->GetXaxis()->SetTitle("cos #theta");
	grCoreNorm->GetXaxis()->SetTitleSize(0.085);
	grCoreNorm->GetXaxis()->SetTitleOffset(0.85);
	grCoreNorm->GetXaxis()->CenterTitle();
	grCoreNorm->GetXaxis()->SetLabelSize(0.065);
	grCoreNorm->GetXaxis()->SetNdivisions(5);
	grCoreNorm->GetXaxis()->SetLimits(-1.0,1.0);
	grCoreNorm->SetMarkerColor(1);
	grCoreNorm->SetMarkerStyle(21);
	grCoreNorm->SetMarkerSize(1.5);
  //grCoreNorm->Draw("AP");

  
  grCoreNorm->Fit("corr420");
  cout << "chisq/ndf 4->2->0: " << corr420->GetChisquare()/(numBins-3) << endl;
  grCoreNorm->Fit("corr520");
  cout << "chisq/ndf 5->2->0: " << corr520->GetChisquare()/(numBins-3) << endl;

  acList->Add(grCoreNorm);
  
  //FILE *out = fopen(outfile, "w");
  //printf("File %s opened\n", outfile);

  cout << "Writing histograms to " << outfile << endl;
  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  acList->Write();

  myfile->Write();
  myfile->Close();


  fclose(inp);
}
int main(int argc, char **argv)
{

  AngCorr *mysort = new AngCorr();

  const char *sfile;
  const char *outfile;
  printf("Starting AngCorrSMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online TIP+TIGRESS data" << endl;
    cout << "Arguments: AngCorrSMOL smol_file eLow1 eHigh1 eLow2 eHigh2 eBGLow1 eBGHigh1 eBGLow2 eBGHigh2 output_file" << endl;
    return 0;
  }else if(argc == 10){
    sfile = argv[1];
    eGamma1[0] = atoi(argv[2]);
    eGamma1[1] = atoi(argv[3]);
    eGamma2[0] = atoi(argv[4]);
    eGamma2[1] = atoi(argv[5]);
    eBG1[0] = atoi(argv[6]);
    eBG1[1] = atoi(argv[7]);
    eBG2[0] = atoi(argv[8]);
    eBG2[1] = atoi(argv[9]);
    outfile = "angcorr.root";
    printf("SMOL file: %s\nOutput file: %s\n", sfile, outfile); 
  }else if(argc == 11){
    sfile = argv[1];
    eGamma1[0] = atoi(argv[2]);
    eGamma1[1] = atoi(argv[3]);
    eGamma2[0] = atoi(argv[4]);
    eGamma2[1] = atoi(argv[5]);
    eBG1[0] = atoi(argv[6]);
    eBG1[1] = atoi(argv[7]);
    eBG2[0] = atoi(argv[8]);
    eBG2[1] = atoi(argv[9]);
    outfile = argv[10];
    printf("SMOL file: %s\nOutput file: %s\n", sfile, outfile);
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: AngCorrSMOL smol_file eLow1 eHigh1 eLow2 eHigh2 eBGLow1 eBGHigh1 eBGLow2 eBGHigh2 output_file\n");
    return 0;
  }

  if(eGamma1[0] > eGamma1[1]){
    //swap values
    int swapVal = eGamma1[1];
    eGamma1[1] = eGamma1[0];
    eGamma1[0] = swapVal;
  }
  if(eGamma2[0] > eGamma2[1]){
    //swap values
    int swapVal = eGamma2[1];
    eGamma2[1] = eGamma2[0];
    eGamma2[0] = swapVal;
  }
  if(eBG1[0] > eBG1[1]){
    //swap values
    int swapVal = eBG1[1];
    eBG1[1] = eBG1[0];
    eBG1[0] = swapVal;
  }
  if(eBG2[0] > eBG2[1]){
    //swap values
    int swapVal = eBG2[1];
    eBG2[1] = eBG2[0];
    eBG2[0] = swapVal;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(sfile, outfile);

  return 0;
}
