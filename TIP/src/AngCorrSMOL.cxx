//use the Makefile!

#define AngCorrS_cxx
#include "common.h"
#include "AngCorrSMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];
int eGamma1[2], eGamma2[2];
double xCore[NUM_HIST_BINS], yCoreRaw[NUM_HIST_BINS], yCoreEvtMix[NUM_HIST_BINS];
double eCoreRaw[NUM_HIST_BINS], eCoreEvtMix[NUM_HIST_BINS];
double yCoreNorm[NUM_HIST_BINS], eCoreNorm[NUM_HIST_BINS];
double xSeg[NUM_HIST_BINS], ySegRaw[NUM_HIST_BINS], ySegEvtMix[NUM_HIST_BINS];
double eSegRaw[NUM_HIST_BINS], eSegEvtMix[NUM_HIST_BINS];
double ySegNorm[NUM_HIST_BINS], eSegNorm[NUM_HIST_BINS];


void AngCorr::SortData(char const *sfile, char const *outfile, const int doppCorr)
{
  Initialise();

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

  float hit1E, hit2E;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(doppCorr == 0){
        hit1E = sortedEvt.noABHit[tigHitInd].energy;
      }else{
        hit1E = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
      }

      if((hit1E >= eGamma1[0])&&(hit1E <= eGamma1[1])){
        for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
          if(sortedEvt.noABHit[tigHitInd].core/4 != sortedEvt.noABHit[tigHitInd2].core/4){ //omit same clover events, to remove cross-talk effects
            if(doppCorr == 0){
              hit2E = sortedEvt.noABHit[tigHitInd2].energy;
            }else{
              hit2E = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
            }
            if((hit2E >= eGamma2[0])&&(hit2E <= eGamma2[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrRaw->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreRaw->Fill(cos(coreangle));
            }
          }
        }
        memcpy(&pastHit[pastEvtInd],&sortedEvt.noABHit[tigHitInd],sizeof(tig_hit));
        pastHitEvt[pastEvtInd] = jentry;
        pastEvtInd++;
        if(pastEvtInd >= NUM_PAST_HIT){
          pastEvtInd = 0; //roll over
        }
      }else if((hit1E >= eGamma2[0])&&(hit1E <= eGamma2[1])){
        for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
          if(sortedEvt.noABHit[tigHitInd].core/4 != sortedEvt.noABHit[tigHitInd2].core/4){ //omit same clover events, to remove cross-talk effects
            if(doppCorr == 0){
              hit2E = sortedEvt.noABHit[tigHitInd2].energy;
            }else{
              hit2E = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitInd2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
            }
            if((hit2E >= eGamma1[0])&&(hit2E <= eGamma1[1])){
              //cascade identified
              TVector3 vec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg);
              TVector3 vec2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,sortedEvt.noABHit[tigHitInd2].seg);
              Double_t angle = vec1.Angle(vec2);
              angCorrRaw->Fill(cos(angle));
              TVector3 corevec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
              TVector3 corevec2 = getTigVector(sortedEvt.noABHit[tigHitInd2].core,0);
              Double_t coreangle = corevec1.Angle(corevec2);
              angCorrCoreRaw->Fill(cos(coreangle));
            }
          }
        }
        memcpy(&pastHit[pastEvtInd],&sortedEvt.noABHit[tigHitInd],sizeof(tig_hit));
        pastHitEvt[pastEvtInd] = jentry;
        pastEvtInd++;
        if(pastEvtInd >= NUM_PAST_HIT){
          pastEvtInd = 0; //roll over
        }
      }

      //handle event mixing
      if((hit1E >= eGamma1[0])&&(hit1E <= eGamma1[1])){
        for(int tigHitInd2 = 0; tigHitInd2 < NUM_PAST_HIT; tigHitInd2++){
          if((pastHitEvt[tigHitInd2] < jentry)&&(pastHitEvt[tigHitInd2] > (jentry-1000))){
            if(sortedEvt.noABHit[tigHitInd].core/4 != pastHit[tigHitInd2].core/4){ //omit same clover events, to remove cross-talk effects
              //hit is from a different event
              if(doppCorr == 0){
                hit2E = pastHit[tigHitInd2].energy;
              }else{
                hit2E = getEDoppFusEvapDirect(&pastHit[tigHitInd2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
              }
              if((hit2E >= eGamma2[0])&&(hit2E <= eGamma2[1])){
                //event mixed cascade identified
                TVector3 vec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg);
                TVector3 vec2 = getTigVector(pastHit[tigHitInd2].core,pastHit[tigHitInd2].seg);
                Double_t angle = vec1.Angle(vec2);
                angCorrEvtMix->Fill(cos(angle));
                TVector3 corevec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
                TVector3 corevec2 = getTigVector(pastHit[tigHitInd2].core,0);
                Double_t coreangle = corevec1.Angle(corevec2);
                angCorrCoreEvtMix->Fill(cos(coreangle));
                evtMixDist->Fill(jentry - pastHitEvt[tigHitInd2]);
              }
            }
          }
        }
      }else if((hit1E >= eGamma2[0])&&(hit1E <= eGamma2[1])){
        for(int tigHitInd2 = 0; tigHitInd2 < NUM_PAST_HIT; tigHitInd2++){
          if((pastHitEvt[tigHitInd2] < jentry)&&(pastHitEvt[tigHitInd2] > (jentry-1000))){
            if(sortedEvt.noABHit[tigHitInd].core/4 != pastHit[tigHitInd2].core/4){ //omit same clover events, to remove cross-talk effects
              //hit is from a different event
              if(doppCorr == 0){
                hit2E = pastHit[tigHitInd2].energy;
              }else{
                hit2E = getEDoppFusEvapDirect(&pastHit[tigHitInd2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
              }
              if((hit2E >= eGamma1[0])&&(hit2E <= eGamma1[1])){
                //event mixed cascade identified
                TVector3 vec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,sortedEvt.noABHit[tigHitInd].seg);
                TVector3 vec2 = getTigVector(pastHit[tigHitInd2].core,pastHit[tigHitInd2].seg);
                Double_t angle = vec1.Angle(vec2);
                angCorrEvtMix->Fill(cos(angle));
                TVector3 corevec1 = getTigVector(sortedEvt.noABHit[tigHitInd].core,0);
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

  //dump histograms to arrays (core)
  int bin=0, numBins=0;
  xCore[bin]=0.;
  yCoreRaw[bin]=0.;
  eCoreRaw[bin]=0.;
  yCoreEvtMix[bin]=0.;
  eCoreEvtMix[bin]=0.;
  cout << "Bins: [ " << angCorrCoreRaw->GetNbinsX() << " " << angCorrCoreEvtMix->GetNbinsX() << " ]." << endl;
  for(int i=0; i<angCorrCoreRaw->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrCoreRaw->GetBinContent(i) > 0.){
        //cout << angCorrCoreRaw->GetBinCenter(i) << " " << angCorrCoreRaw->GetBinContent(i) << endl;
        xCore[bin]=angCorrCoreRaw->GetBinCenter(i);
        yCoreRaw[bin]=angCorrCoreRaw->GetBinContent(i);
        eCoreRaw[bin] = sqrt(yCoreRaw[bin]);
        bin++;
        xCore[bin]=0.;
        yCoreRaw[bin]=0.;
        eCoreRaw[bin]=0.;
      }
    }
  }
  numBins=bin;
  bin=0;
  for(int i=0; i<angCorrCoreEvtMix->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrCoreEvtMix->GetBinContent(i) > 0.){
        if(angCorrCoreEvtMix->GetBinCenter(i) != angCorrCoreRaw->GetBinCenter(i)){
          printf("x-value mismatch between raw and event mixed at bin %i: %f %f\n",i,xCore[bin],angCorrCoreEvtMix->GetBinCenter(i));
        }
        yCoreEvtMix[bin]+=angCorrCoreEvtMix->GetBinContent(i);
        eCoreEvtMix[bin] = sqrt(yCoreEvtMix[bin]);
        bin++;
        yCoreEvtMix[bin]=0.;
        eCoreEvtMix[bin]=0.;
      }
    }
  }

  //generate normalized data
  double avg = 0.;
  for(int i=0; i<numBins-1; i++){
    yCoreNorm[i] = (yCoreRaw[i])/yCoreEvtMix[i];
    eCoreNorm[i] = sqrt(pow(eCoreRaw[i]/yCoreRaw[i],2) + pow(eCoreEvtMix[i]/yCoreEvtMix[i],2))*yCoreNorm[i];
    avg += yCoreNorm[i];
  }
  avg /= (numBins-1);
  for(int i=0; i<numBins-1; i++){
    yCoreNorm[i] /= avg;
    eCoreNorm[i] /= avg;
  }

  //printf normalized data
  cout << "cos(theta)    W(theta)    err" << endl;
  for(int i=0; i<numBins-1; i++){
    if(yCoreNorm[i] > 0.){
      cout << xCore[i] << " " << yCoreNorm[i] << " " << eCoreNorm[i] << endl;
    }
  }

  //plot data with error bars 
	TGraphErrors *grCoreNorm = new TGraphErrors(numBins-1,xCore,yCoreNorm,nullptr,eCoreNorm);
	grCoreNorm->SetTitle("Normalized angular correlation");
	grCoreNorm->GetYaxis()->SetTitle("W(#theta)");
	grCoreNorm->GetYaxis()->SetTitleSize(0.085);
	grCoreNorm->GetYaxis()->SetTitleOffset(0.5);
	grCoreNorm->GetYaxis()->CenterTitle();
	grCoreNorm->GetYaxis()->SetLabelSize(0.065);
	grCoreNorm->GetYaxis()->SetNdivisions(6);
	grCoreNorm->GetYaxis()->SetDecimals();
	grCoreNorm->GetYaxis()->SetRangeUser(0.0,2.0);
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

  
  grCoreNorm->Fit("corrfit");
  cout << "chisq/ndf: " << corrfit->GetChisquare()/(numBins-1) << endl;
  grCoreNorm->Fit("corrfitzeroa4");
  cout << "chisq/ndf (a4 = 0): " << corrfitzeroa4->GetChisquare()/(numBins-1) << endl;

  acList->Add(grCoreNorm);

  //dump histograms to arrays (segments)
  bin=0, numBins=0;
  xSeg[bin]=0.;
  ySegRaw[bin]=0.;
  eSegRaw[bin]=0.;
  ySegEvtMix[bin]=0.;
  eSegEvtMix[bin]=0.;
  cout << "Bins: [ " << angCorrRaw->GetNbinsX() << " " << angCorrEvtMix->GetNbinsX() << " ]." << endl;
  for(int i=0; i<angCorrRaw->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrRaw->GetBinContent(i) > 0.){
        //cout << angCorrRaw->GetBinCenter(i) << " " << angCorrRaw->GetBinContent(i) << endl;
        xSeg[bin]=angCorrRaw->GetBinCenter(i);
        ySegRaw[bin]=angCorrRaw->GetBinContent(i);
        eSegRaw[bin] = sqrt(ySegRaw[bin]);
        bin++;
        xSeg[bin]=0.;
        ySegRaw[bin]=0.;
        eSegRaw[bin]=0.;
      }
    }
  }
  numBins=bin;
  bin=0;
  for(int i=0; i<angCorrEvtMix->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(angCorrEvtMix->GetBinContent(i) > 0.){
        if(angCorrEvtMix->GetBinCenter(i) != angCorrRaw->GetBinCenter(i)){
          printf("x-value mismatch between raw and event mixed at bin %i: %f %f\n",i,xSeg[bin],angCorrEvtMix->GetBinCenter(i));
        }
        ySegEvtMix[bin]+=angCorrEvtMix->GetBinContent(i);
        eSegEvtMix[bin] = sqrt(ySegEvtMix[bin]);
        bin++;
        ySegEvtMix[bin]=0.;
        eSegEvtMix[bin]=0.;
      }
    }
  }

  //generate normalized data
  avg = 0.;
  for(int i=0; i<numBins-1; i++){
    ySegNorm[i] = (ySegRaw[i])/ySegEvtMix[i];
    eSegNorm[i] = sqrt(pow(eSegRaw[i]/ySegRaw[i],2) + pow(eSegEvtMix[i]/ySegEvtMix[i],2))*ySegNorm[i];
    avg += ySegNorm[i];
  }
  avg /= (numBins-1);
  for(int i=0; i<numBins-1; i++){
    ySegNorm[i] /= avg;
    eSegNorm[i] /= avg;
  }

  //printf normalized data
  cout << "cos(theta)    W(theta)    err" << endl;
  for(int i=0; i<numBins-1; i++){
    if(ySegNorm[i] > 0.){
      cout << xSeg[i] << " " << ySegNorm[i] << " " << eSegNorm[i] << endl;
    }
  }

  //plot data with error bars 
	TGraphErrors *grSegNorm = new TGraphErrors(numBins-1,xSeg,ySegNorm,nullptr,eSegNorm);
	grSegNorm->SetTitle("Normalized angular correlation");
	grSegNorm->GetYaxis()->SetTitle("W(#theta)");
	grSegNorm->GetYaxis()->SetTitleSize(0.085);
	grSegNorm->GetYaxis()->SetTitleOffset(0.5);
	grSegNorm->GetYaxis()->CenterTitle();
	grSegNorm->GetYaxis()->SetLabelSize(0.065);
	grSegNorm->GetYaxis()->SetNdivisions(6);
	grSegNorm->GetYaxis()->SetDecimals();
	grSegNorm->GetYaxis()->SetRangeUser(0.0,2.0);
	grSegNorm->GetXaxis()->SetTitle("cos #theta");
	grSegNorm->GetXaxis()->SetTitleSize(0.085);
	grSegNorm->GetXaxis()->SetTitleOffset(0.85);
	grSegNorm->GetXaxis()->CenterTitle();
	grSegNorm->GetXaxis()->SetLabelSize(0.065);
	grSegNorm->GetXaxis()->SetNdivisions(5);
	grSegNorm->GetXaxis()->SetLimits(-1.0,1.0);
	grSegNorm->SetMarkerColor(1);
	grSegNorm->SetMarkerStyle(21);
	grSegNorm->SetMarkerSize(1.5);
  //grSegNorm->Draw("AP");

  
  grSegNorm->Fit("corrfit");
  cout << "chisq/ndf: " << corrfit->GetChisquare()/(numBins-1) << endl;
  grSegNorm->Fit("corrfitzeroa4");
  cout << "chisq/ndf (a4 = 0): " << corrfitzeroa4->GetChisquare()/(numBins-1) << endl;

  acList->Add(grSegNorm);
  
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
  int doppCorr = 0;
  printf("Starting AngCorrSMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online TIP+TIGRESS data" << endl;
    cout << "Arguments: AngCorrSMOL smol_file doppCorr eLow1 eHigh1 eLow2 eHigh2 utput_file" << endl;
    cout << "doppCorr = 1 to apply Doppler correction" << endl;
    return 0;
  }else if(argc == 7){
    sfile = argv[1];
    doppCorr = atoi(argv[2]);
    eGamma1[0] = atoi(argv[3]);
    eGamma1[1] = atoi(argv[4]);
    eGamma2[0] = atoi(argv[5]);
    eGamma2[1] = atoi(argv[6]);
    outfile = "angcorr.root";
  }else if(argc == 8){
    sfile = argv[1];
    doppCorr = atoi(argv[2]);
    eGamma1[0] = atoi(argv[3]);
    eGamma1[1] = atoi(argv[4]);
    eGamma2[0] = atoi(argv[5]);
    eGamma2[1] = atoi(argv[6]);
    outfile = argv[7];
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: AngCorrSMOL smol_file eLow1 eHigh1 eLow2 eHigh2 output_file\n");
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

  gates = new PIDGates;
  printf("SMOL file: %s\nOutput file: %s\n", sfile, outfile);
  printf("Energy gates: [%i %i] keV, [%i %i] keV\n",eGamma1[0],eGamma1[1],eGamma2[0],eGamma2[1]);
  if(doppCorr == 0){
    cout << "Not performing Doppler correction." << endl;
  }else{
    cout << "Performing Doppler correction." << endl;
  }
  cout << endl;

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(sfile, outfile, doppCorr);

  return 0;
}
