//use the Makefile!

#define AngCorrS_cxx
#include "common.cxx"
#include "AngCorrSMOL.h"

using namespace std;

sorted_evt sortedEvt, prevEvtBuf[MAX_EVT_MIX_OFFSET];

//angular correlation energy gates
Double_t gamma1EGate[2], gamma2EGate[2], gamma2BGEGate[2];

//timing windows
const Double_t hpgehpgeTGate[2] = {-30, 30}; // narrow HPGe - HPGe timing window (ns)
const Double_t hpgehpgeTSGate[2] = {0, 6}; // narrow HPGe - HPGe timing window (timestamp units)
const Double_t hpgehpgeTRandGate[2] = {1300, 2000}; // time-random HPGe - HPGe timing window (ns)

void WriteData(const char* outName){

  cout << "Writing histograms to " << outName << endl;

  TFile *myfile = new TFile(outName, "RECREATE");
  myfile->cd();
  TDirectory *angcorrdir = myfile->mkdir("AngularCorrelation");
  angcorrdir->cd();
  angcorrList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();

}


void SortData(const char *sfile, const uint8_t forwardPos)
{

  FILE *inp = fopen(sfile, "rb");
  if(inp==NULL){
    printf("ERROR: couldn't open file: %s\n",sfile);
    exit(-1);
  }
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  uint64_t pileupCtrs[16];
  fread(&sentries,sizeof(uint64_t),1,inp);
  uint64_t smolVersion = (uint64_t)(sentries >> 48);
  if(smolVersion > 0){
    fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
    //printf("\nNumber of hits of each pileup type:\n");
    uint64_t totalHits = 0;
    for(uint8_t i=0; i<16; i++){
      //printf("Pileup type %2u: %Lu\n",i,pileupCtrs[i]);
      totalHits += pileupCtrs[i];
    }
    printf("Total hits:     %Lu\n",totalHits);
    long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
    printf("Fraction of hits with no pileup: %Lf\n",frac);
  }
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  
  memset(prevEvtBuf,0,sizeof(prevEvtBuf));
  int64_t prevEvtBufPos = 0; //buffer position that we're currently at
  int64_t numPrevEvtBufEntries = 0;
  
  uint8_t footerVal;

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
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        if((sortedEvt.noABHit[noABHitInd].energy >= gamma1EGate[0])&&(sortedEvt.noABHit[noABHitInd].energy <= gamma1EGate[1])){
          for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
            if(sortedEvt.noABHit[noABHitInd2].energy > MIN_HPGE_EAB){
              Double_t tDiff = (noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd));
              if((tDiff >= COINC_TIMING_GATE_MIN)&&(tDiff <= COINC_TIMING_GATE_MAX)){
                if((sortedEvt.noABHit[noABHitInd2].energy >= gamma2EGate[0])&&(sortedEvt.noABHit[noABHitInd2].energy <= gamma2EGate[1])){
                  Double_t theta = getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos));
                  angCorrRaw->Fill(cos(theta));
                }else if((sortedEvt.noABHit[noABHitInd2].energy >= gamma2BGEGate[0])&&(sortedEvt.noABHit[noABHitInd2].energy <= gamma2BGEGate[1])){
                  Double_t theta = getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(sortedEvt.noABHit[noABHitInd2].core & 63U,forwardPos));
                  angCorrCompton->Fill(cos(theta));
                }
              }
            }
          }
        }
      }
    }

    //do event mixing if data is available in the previous event buffer
    if(numPrevEvtBufEntries >= MIN_EVT_MIX_OFFSET){
      int64_t evtToMixInd;
      for(int64_t i=MIN_EVT_MIX_OFFSET; i<numPrevEvtBufEntries; i++){
        if(prevEvtBufPos >= i){
          evtToMixInd = prevEvtBufPos - i;
        }else{
          evtToMixInd = (numPrevEvtBufEntries)-(i-prevEvtBufPos);
        }
        if(((prevEvtBufPos - evtToMixInd) >= 0)&&((prevEvtBufPos - evtToMixInd) < MIN_EVT_MIX_OFFSET)){
          printf("i: %i\n",i);
          printf(" numPrevEvtBufEntries: %i\n",numPrevEvtBufEntries);
          printf(" prevEvtBufPos: %i\n",prevEvtBufPos);
          printf(" evtToMixInd: %i\n",evtToMixInd);
          printf("offset: %i\n",prevEvtBufPos - evtToMixInd);
        }

        for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
          if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
            if((sortedEvt.noABHit[noABHitInd].energy >= gamma1EGate[0])&&(sortedEvt.noABHit[noABHitInd].energy <= gamma1EGate[1])){
              for(int noABHitInd2 = 0; noABHitInd2 < prevEvtBuf[evtToMixInd].header.numNoABHits; noABHitInd2++){
                if(prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].energy > MIN_HPGE_EAB){
                  //no coincidence timing condition since we are looking at different events
                  if((prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].energy >= gamma2EGate[0])&&(prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].energy <= gamma2EGate[1])){
                    Double_t theta = getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].core & 63U,forwardPos));
                    angCorrEvtMix->Fill(cos(theta));
                  }else if((prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].energy >= gamma2BGEGate[0])&&(prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].energy <= gamma2BGEGate[1])){
                    Double_t theta = getGRIFFINVector(sortedEvt.noABHit[noABHitInd].core & 63U,forwardPos).Angle(getGRIFFINVector(prevEvtBuf[evtToMixInd].noABHit[noABHitInd2].core & 63U,forwardPos));
                    angCorrEvtMixCompton->Fill(cos(theta));
                  }
                }
              }
            }
          }
        }
      }
    }


    //copy event data to previous event buffer, to allow for event mixing
    memcpy(&prevEvtBuf[prevEvtBufPos],&sortedEvt,sizeof(sorted_evt));
    prevEvtBufPos++;
    if(prevEvtBufPos >= MAX_EVT_MIX_OFFSET){
      prevEvtBufPos = 0; //circle back
    }
    if(numPrevEvtBufEntries < MAX_EVT_MIX_OFFSET){
      numPrevEvtBufEntries++; //buffer is being filled
    }

    if (jentry % 9713 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;

  fclose(inp);
}
int main(int argc, char **argv)
{

  const char *sfile;
  const char *outfile;
  char outTxtFilename[256], tmpStr[256];
  uint8_t forwardPos = 0;
  printf("Starting AngCorrSMOL\n");

  if(argc < 8){
    cout << "Code sorts angular correlation ROOT histograms for the specified cascade." << endl;
    cout << "Arguments: AngCorrSMOL smol_file gamma_energy_1 gamma_energy_2 BG_energy gate_width forward_pos output_file" << endl;
    cout << "  *smol_file* can be a single SMOL tree (extension .smol), or a list of SMOL trees (extension .list, one filepath per line)." << endl;
    cout << "  *gamma_energy_1* and *gamma_energy_2* are energies of the 2 gammas in the cascade, in keV." << endl;
    cout << "  *BG_energy* is an energy in keV correponding to Compton background near the 2nd gamma, for background subtraction." << endl;
    cout << "  *gate_width* is the range in keV around each specified energy, and should be approximately the peak widths." << endl;
    cout << "  *forward_pos* should 1 if GRIFFIN is at 110 mm, 0 if at 145 mm." << endl;
    return 0;
  }else if(argc > 8){
    printf("ERROR: too many arguments!\nArguments: AngCorrSMOL smol_file gamma_energy_1 gamma_energy_2 BG_energy gate_width output_file\n");
    return 0;
  }

  sfile = argv[1];
  outfile = argv[7];
  Double_t gamma1E = atof(argv[2]);
  Double_t gamma2E = atof(argv[3]);
  Double_t gamma2BGE = atof(argv[4]);
  Double_t gateWidth = atof(argv[5]);
  forwardPos = (uint8_t)atoi(argv[6]);
  if((gamma1E < 0)||(gamma2E < 0)||(gamma2BGE < 0)||(gateWidth < 0)){
    cout << "ERROR: gamma energies and gate widths must all be positive." << endl;
    return 0;
  }
  gamma1EGate[0] = gamma1E - gateWidth/2.0;
  gamma1EGate[1] = gamma1E + gateWidth/2.0;
  if(gamma1EGate[0] < 0){
    gamma1EGate[0] = 0.0;
  }
  gamma2EGate[0] = gamma2E - gateWidth/2.0;
  gamma2EGate[1] = gamma2E + gateWidth/2.0;
  if(gamma2EGate[0] < 0){
    gamma2EGate[0] = 0.0;
  }
  gamma2BGEGate[0] = gamma2BGE - gateWidth/2.0;
  gamma2BGEGate[1] = gamma2BGE + gateWidth/2.0;
  if(gamma2BGEGate[0] < 0){
    gamma2BGEGate[0] = 0.0;
  }

  printf("Energy gates:\n");
  printf("  Gamma 1: [%.2f %.2f] keV\n", gamma1EGate[0], gamma1EGate[1]);
  printf("  Gamma 2: [%.2f %.2f] keV\n", gamma2EGate[0], gamma2EGate[1]);
  printf("  Gamma 2 background: [%.2f %.2f] keV\n", gamma2BGEGate[0], gamma2BGEGate[1]);
  if(forwardPos){
    printf("GRIFFIN at 110 mm.\n");
  }else{
    printf("GRIFFIN at 145 mm.\n");
  }

  theApp=new TApplication("App", &argc, argv);
  InitialiseHists();

  strncpy(tmpStr,outfile,256);
  char *tok = strtok(tmpStr,".");
  snprintf(outTxtFilename,256,"%s.txt",tok);
  printf("Will save angular distribution plaintext data to: %s\n",outTxtFilename);

  const char *dot = strrchr(sfile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get SMOL tree or list file extension." << endl;
    return 0;
  }

  uint64_t numSepEvts = 0U;
  if(strcmp(dot + 1, "smol") == 0){
    printf("SMOL tree: %s\nOutput file: %s\n", sfile, outfile);
    SortData(sfile,forwardPos);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("SMOL tree list: %s\nOutput file: %s\n", sfile, outfile);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(sfile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << sfile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          SortData(str,forwardPos);
        }
      }
    }
  }else{
    printf("ERROR: improper file extension for SMOL tree or list (should be .smol or .list).\n");
    return 0;
  }

  FILE *outTxt = fopen(outTxtFilename,"w");


  //make event mixing corrected distribution
  int numAngCorrDataPts = 0;
  int numAngCorrDataPtsCompton = 0;
  Double_t numerErr, denomErr;
  Double_t angCorrSum = 0.0;
  Double_t angCorrSumCompton = 0.0;
  for(int i=0; i<=angCorrRaw->GetNbinsX(); i++){
    if((angCorrRaw->GetBinContent(i) > 0)&&(angCorrEvtMix->GetBinContent(i) > 0)){
      angCorrEvtMixCorr->SetBinContent(i,angCorrRaw->GetBinContent(i)/angCorrEvtMix->GetBinContent(i));
      numerErr = angCorrRaw->GetBinError(i);
      denomErr = angCorrEvtMix->GetBinError(i);
      angCorrEvtMixCorr->SetBinError(i,(angCorrEvtMixCorr->GetBinContent(i))*sqrt(pow(numerErr/angCorrRaw->GetBinContent(i),2.0) + pow(denomErr/angCorrEvtMix->GetBinContent(i),2.0)));
      angCorrSum += angCorrEvtMixCorr->GetBinContent(i);
      numAngCorrDataPts++;
    }
    if(((angCorrRaw->GetBinContent(i) - angCorrCompton->GetBinContent(i)) > 0)&&((angCorrEvtMix->GetBinContent(i) - angCorrEvtMixCompton->GetBinContent(i)) > 0)){
      angCorrEvtMixCorrCompton->SetBinContent(i,(angCorrRaw->GetBinContent(i) - angCorrCompton->GetBinContent(i))/(angCorrEvtMix->GetBinContent(i) - angCorrEvtMixCompton->GetBinContent(i)));
      numerErr = sqrt(pow(angCorrRaw->GetBinError(i),2.0) + pow(angCorrCompton->GetBinError(i),2.0));
      denomErr = sqrt(pow(angCorrEvtMix->GetBinError(i),2.0) + pow(angCorrEvtMixCompton->GetBinError(i),2.0));
      angCorrEvtMixCorrCompton->SetBinError(i,(angCorrEvtMixCorr->GetBinContent(i))*sqrt(pow(numerErr/(angCorrRaw->GetBinContent(i) - angCorrCompton->GetBinContent(i)),2.0) + pow(denomErr/(angCorrEvtMix->GetBinContent(i) - angCorrEvtMixCompton->GetBinContent(i)),2.0)));
      angCorrSumCompton += angCorrEvtMixCorrCompton->GetBinContent(i);
      numAngCorrDataPtsCompton++;
    }
  }
  //normalize and print angular correlation
  if(angCorrSum > 0){
    fprintf(outTxt,"Event mixing corrected angular correlation:\n");
    fprintf(outTxt,"cos(theta) ang_corr_value error\n");
    Double_t scaling = (1.0*numAngCorrDataPts)/angCorrSum;
    for(int i=0; i<=angCorrEvtMixCorr->GetNbinsX(); i++){
      Double_t err = angCorrEvtMixCorr->GetBinError(i);
      angCorrEvtMixCorr->SetBinContent(i,angCorrEvtMixCorr->GetBinContent(i)*scaling);
      angCorrEvtMixCorr->SetBinError(i,err*scaling);
      if(angCorrEvtMixCorr->GetBinContent(i) > 0.01){
        fprintf(outTxt,"%0.5f %0.5f %0.5f\n",angCorrEvtMixCorr->GetBinCenter(i),angCorrEvtMixCorr->GetBinContent(i),angCorrEvtMixCorr->GetBinError(i));
      }
    }
  }else{
    printf("ERROR: angular correlation is negative!\n");
    //return 0;
  }
  if(angCorrSumCompton > 0){
    fprintf(outTxt,"\nEvent mixing and Compton corrected angular correlation:\n");
    fprintf(outTxt,"cos(theta) ang_corr_value error\n");
    Double_t scaling = (1.0*numAngCorrDataPtsCompton)/angCorrSumCompton;
    for(int i=0; i<=angCorrEvtMixCorrCompton->GetNbinsX(); i++){
      Double_t err = angCorrEvtMixCorrCompton->GetBinError(i);
      angCorrEvtMixCorrCompton->SetBinContent(i,angCorrEvtMixCorrCompton->GetBinContent(i)*scaling);
      angCorrEvtMixCorrCompton->SetBinError(i,err*scaling);
      if(angCorrEvtMixCorrCompton->GetBinContent(i) > 0.01){
        fprintf(outTxt,"%0.5f %0.5f %0.5f\n",angCorrEvtMixCorrCompton->GetBinCenter(i),angCorrEvtMixCorrCompton->GetBinContent(i),angCorrEvtMixCorrCompton->GetBinError(i));
      }
    }
  }else{
    printf("ERROR: Compton corrected angular correlation is negative!\n");
    //return 0;
  }

  fclose(outTxt);
  WriteData(outfile);

  return 0;
}
