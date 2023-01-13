//use the Makefile!

#define AngDistS_cxx
#include "common.h"
#include "AngDistSMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];
int eGamma[2], eGammaBG[2], eGammaSource[2];
double x[NUM_HIST_BINS], yRaw[NUM_HIST_BINS], yBG[NUM_HIST_BINS], ySource[NUM_HIST_BINS];
double eRaw[NUM_HIST_BINS], eBG[NUM_HIST_BINS], eSource[NUM_HIST_BINS];
double yNorm[NUM_HIST_BINS], eNorm[NUM_HIST_BINS];


void AngDist::SortData(char const *sfile, char const *sfilesrc, char const *outfile)
{
  Initialise();

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nSorting events...\n");
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

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigABHits; tigHitInd++){
      float addE1 = sortedEvt.tigHit[tigHitInd].energy;

      if(addE1 > MIN_TIG_EAB){
        if((addE1 >= eGamma[0])&&(addE1 <= eGamma[1])){
          //gamma identified
          TVector3 vec = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
          Double_t angle = vec.Theta();
          AngDistRaw->Fill(cos(angle));
          TVector3 corevec = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
          Double_t coreangle = corevec.Theta();
          AngDistCoreRaw->Fill(cos(coreangle));
        }
        if((addE1 >= eGammaBG[0])&&(addE1 <= eGammaBG[1])){
          //gamma identified
          TVector3 vec = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
          Double_t angle = vec.Theta();
          AngDistBG->Fill(cos(angle));
          TVector3 corevec = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
          Double_t coreangle = corevec.Theta();
          AngDistCoreBG->Fill(cos(coreangle));
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;
  fclose(inp);


  inp = fopen(sfilesrc, "rb");
  printf("File %s opened\n", sfilesrc);
  
  sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);

  printf("\nSorting events...\n");
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

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigABHits; tigHitInd++){
      float addE1 = sortedEvt.tigHit[tigHitInd].energy;

      if(addE1 > MIN_TIG_EAB){
        if((addE1 >= eGammaSource[0])&&(addE1 <= eGammaSource[1])){
          //gamma identified
          TVector3 vec = getTigVector(sortedEvt.tigHit[tigHitInd].core,sortedEvt.tigHit[tigHitInd].seg);
          Double_t angle = vec.Theta();
          AngDistSource->Fill(cos(angle));
          TVector3 corevec = getTigVector(sortedEvt.tigHit[tigHitInd].core,0);
          Double_t coreangle = corevec.Theta();
          AngDistCoreSource->Fill(cos(coreangle));
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;
  fclose(inp);


  //dump histograms to arrays
  int bin=0, numBins=0;
  for(int i=0; i<AngDistRaw->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(AngDistRaw->GetBinContent(i) > 0.){
        x[bin]=AngDistRaw->GetBinCenter(i);
        yRaw[bin]=AngDistRaw->GetBinContent(i);
        eRaw[bin]=AngDistRaw->GetBinError(i);
        bin++;
      }
    }
  }
  numBins=bin;
  bin=0;
  for(int i=0; i<AngDistBG->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      if(AngDistRaw->GetBinContent(i) > 0.){
        if(AngDistBG->GetBinCenter(i) != x[bin]){
          printf("x-value mismatch between raw and BG at bin %i: %f %f\n",i,x[bin],AngDistBG->GetBinCenter(i));
          bin++; //try to correct
        }
        yBG[bin]=AngDistBG->GetBinContent(i);
        eBG[bin]=AngDistBG->GetBinError(i);
        bin++;
      }
    }
  }
  if(bin!=numBins){
    printf("Bin number mismatch between raw and BG: %i %i\n",numBins,bin);
  }
  bin=0;
  for(int i=0; i<AngDistSource->GetNbinsX(); i++){
    if(bin<NUM_HIST_BINS){
      //I think some of the segments must be missing in the source data
      //compared to the run data, using a threshold here
      if(AngDistSource->GetBinContent(i) > 20.){
        if(AngDistSource->GetBinCenter(i) != x[bin]){
          printf("x-value mismatch between raw and event mixed at bin %i: %f %f\n",i,x[bin],AngDistSource->GetBinCenter(i));
          //try to correct
          for(int j=bin;j<numBins-1;j++){
            x[j]=x[j+1];
            yRaw[j]=yRaw[j+1];
            eRaw[j]=eRaw[j+1];
            yBG[j]=yBG[j+1];
            eBG[j]=eBG[j+1];
          }
          numBins--;
        }
        ySource[bin]=AngDistSource->GetBinContent(i);
        eSource[bin]=AngDistSource->GetBinError(i);
        bin++;
      }
    }
  }
  if(bin!=numBins){
    printf("Bin number mismatch between raw and BG: %i %i\n",numBins,bin);
  }

  //generate normalized data
  double avg = 0.;
  for(int i=0; i<numBins; i++){
    if(ySource[i] > 20.){
      cout << yRaw[i] << " " << yBG[i] << " " << ySource[i] << endl;
      yNorm[i] = (yRaw[i] - yBG[i])/ySource[i];
      eNorm[i] = ((sqrt(pow(eRaw[i],2) + pow(eBG[i],2)))/(yRaw[i] - yBG[i]));
      eNorm[i] += eSource[i]/ySource[i];
      eNorm[i] *= yNorm[i];
      avg += yNorm[i];
    }
  }
  avg /= numBins;
  for(int i=0; i<numBins; i++){
    yNorm[i] /= avg;
    eNorm[i] /= avg;
    cout << yNorm[i] << endl;
  }

  //plot data with error bars 
	TGraphErrors *grNorm = new TGraphErrors(numBins,x,yNorm,nullptr,eNorm);
	grNorm->SetTitle("Normalized angular distribution");
	grNorm->GetYaxis()->SetTitle("W(#theta)");
	grNorm->GetYaxis()->SetTitleSize(0.085);
	grNorm->GetYaxis()->SetTitleOffset(0.5);
	grNorm->GetYaxis()->CenterTitle();
	grNorm->GetYaxis()->SetLabelSize(0.065);
	grNorm->GetYaxis()->SetNdivisions(6);
	grNorm->GetYaxis()->SetDecimals();
	grNorm->GetYaxis()->SetRangeUser(0.4,1.6);
	grNorm->GetXaxis()->SetTitle("cos #theta");
	grNorm->GetXaxis()->SetTitleSize(0.085);
	grNorm->GetXaxis()->SetTitleOffset(0.85);
	grNorm->GetXaxis()->CenterTitle();
	grNorm->GetXaxis()->SetLabelSize(0.065);
	grNorm->GetXaxis()->SetNdivisions(5);
	grNorm->GetXaxis()->SetLimits(-1.0,1.0);
	grNorm->SetMarkerColor(1);
	grNorm->SetMarkerStyle(21);
	grNorm->SetMarkerSize(1.5);
  //grNorm->Draw("AP");

  
  /*grNorm->Fit("corr420");
  cout << "chisq/ndf 4->2->0: " << corr420->GetChisquare()/(numBins-(5/contractFactor)) << endl;
  grNorm->Fit("corr520");
  cout << "chisq/ndf 5->2->0: " << corr520->GetChisquare()/(numBins-(5/contractFactor)) << endl;*/

  adList->Add(grNorm);
  
  //FILE *out = fopen(outfile, "w");
  //printf("File %s opened\n", outfile);

  cout << "Writing histograms to " << outfile << endl;
  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  adList->Write();

  myfile->Write();
  myfile->Close();


  
}
int main(int argc, char **argv)
{

  AngDist *mysort = new AngDist();

  const char *sfile, *sfilesrc;
  const char *outfile;
  printf("Starting AngDistSMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online TIP+TIGRESS data" << endl;
    cout << "Arguments: AngDistSMOL smol_file eLow eHigh eBGLow eBGHigh smol_file_src eSourceLow eSourceHigh output_file" << endl;
    return 0;
  }else if(argc == 9){
    sfile = argv[1];
    eGamma[0] = atoi(argv[2]);
    eGamma[1] = atoi(argv[3]);
    eGammaBG[0] = atoi(argv[4]);
    eGammaBG[1] = atoi(argv[5]);
    sfilesrc = argv[6];
    eGammaSource[0] = atoi(argv[7]);
    eGammaSource[1] = atoi(argv[8]);
    outfile = "AngDist.root";
    printf("SMOL file (signal): %s\nSMOL file (source): %s\nOutput file: %s\n", sfile, sfilesrc, outfile); 
  }else if(argc == 10){
    sfile = argv[1];
    eGamma[0] = atoi(argv[2]);
    eGamma[1] = atoi(argv[3]);
    eGammaBG[0] = atoi(argv[4]);
    eGammaBG[1] = atoi(argv[5]);
    sfilesrc = argv[6];
    eGammaSource[0] = atoi(argv[7]);
    eGammaSource[1] = atoi(argv[8]);
    outfile = argv[9];
    printf("SMOL file (signal): %s\nSMOL file (source): %s\nOutput file: %s\n", sfile, sfilesrc, outfile); 
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: AngDistSMOL smol_file eLow eHigh eBGLow eBGHigh smol_file_src eSourceLow eSourceHigh output_file\n");
    return 0;
  }

  if(eGamma[0] > eGamma[1]){
    //swap values
    int swapVal = eGamma[1];
    eGamma[1] = eGamma[0];
    eGamma[0] = swapVal;
  }
  if(eGammaBG[0] > eGammaBG[1]){
    //swap values
    int swapVal = eGammaBG[1];
    eGammaBG[1] = eGammaBG[0];
    eGammaBG[0] = swapVal;
  }

  if(eGammaSource[0] > eGammaSource[1]){
    //swap values
    int swapVal = eGammaSource[1];
    eGammaSource[1] = eGammaSource[0];
    eGammaSource[0] = swapVal;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(sfile, sfilesrc, outfile);

  return 0;
}
