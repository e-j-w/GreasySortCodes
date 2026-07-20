//g++ seg_check_writewaveforms.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -o SegCheckWrite
//works with GRSISort commit 213cdaf88ee554cb104385a0d2a9b03da91a1ae8
//and GRSIData commit f5b1c9b3eb2189289c78df6f3ad5a44dc9ca6a39

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TPulseAnalyzer.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

TApplication *theApp;
char histName[256];

void seg_check(char const *infile, char const *calfile, char const *outfile)
{

  FILE *out;
  if((out = fopen(outfile, "w")) == NULL){ //open the file
    cout << "ERROR: Cannot open the output file: " << outfile << endl;
    return;
  }

  TFile *inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()){
    printf("Opening file failed, aborting\n");
    return;
  }
  TChain *AnalysisTree = (TChain *)inputfile->Get("AnalysisTree");
  printf("%i tree files, details:\n", AnalysisTree->GetNtrees());
  AnalysisTree->ls();
  TTree *tree = (TTree *)AnalysisTree->GetTree();
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TTigress *tigr = nullptr;
  if (AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigr);
  }else{
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer!" << endl;
  }

  TTigressHit *tigress_hit;
  int waveform_counter = 0;
  int segwaveform_counter[8];
  for (int i = 0; i < 8; i++)
    segwaveform_counter[i] = 0;
  const std::vector<Short_t> *wf;
  const std::vector<Short_t> *segwf;
  int event_counter = 0;
  double baseline;
  double max;
  double amplitude;
  double tmpTime, tmpTimeSeg, tmpCfd, tmpCfdSeg, tmpCfdT0, tmpCfdSegT0;
  printf("Begin sort\n");
  for (int jentry = 0; jentry < 1000; jentry++)
  {
    tree->GetEntry(jentry);
    for (int one = 0; one < tigr->GetMultiplicity(); one++)
    {

      tigress_hit = tigr->GetTigressHit(one);
      tigress_hit->SetWavefit();
      wf = tigress_hit->GetWaveform();
      if(wf->size() <= 0){ continue; }
      /*baseline = 0;
      for (unsigned int i = 0; i < wf->size(); i++)
      {
        if (i < 20)
          baseline += wf->at(i);
        if (i == 0)
          max = wf->at(i);
        if (max < wf->at(i))
          max = wf->at(i);
      }
      amplitude = (abs(max - baseline * 0.05));
      printf("\nCore\n----\n");
      printf("FV waveform size: %i\nAmplitude: %f",wf->size(),amplitude);*/
      if(wf->size() > 0){
        fprintf(out,"FV: %i,",wf->at(0));
      }
      for (unsigned int i = 1; i < wf->size()-1; i++){
        fprintf(out," %i,",wf->at(i));
      }
      fprintf(out," %i\n",wf->at(wf->size()-1));
      //printf("\nSegments\n--------\n");
      if(tigress_hit->GetSegmentMultiplicity()!=8) continue;

      for (int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++)
      {
        const TDetectorHit* segment_hit = &tigress_hit->GetSegmentHit(j);
        segwf = segment_hit->GetWaveform();
        /*baseline = 0;
        for (unsigned int i = 0; i < segwf->size(); i++)
        {
          if (i < 20)
            baseline += segwf->at(i);
          if (i == 0)
            max = segwf->at(i);
          if (max < segwf->at(i))
            max = segwf->at(i);
        }
        amplitude = (abs(max - baseline * 0.05));
        printf("segment %i waveform size: %i, amplitude: %f\n",segment_hit->GetSegment(),segwf->size(),amplitude);*/
        if(segwf->size() > 0){
          fprintf(out,"%i: %i,",segment_hit->GetSegment(),segwf->at(0));
        }
        for (unsigned int i = 1; i < segwf->size()-1; i++){
          fprintf(out," %i,",segwf->at(i));
        }
        fprintf(out," %i\n",segwf->at(segwf->size()-1));
      }
      waveform_counter++;

      fprintf(out,"\n");
      //getc(stdin);

    }
    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "TIGRESS Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete"
           << "\r" << flush;
  }
  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";
  cout << "Histograms written, sorting complete" << endl;
  cout << "Writing histograms to " << outfile << endl;

  fclose(out);

}

int main(int argc, char **argv)
{

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0)
  {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  // Input-chain-file, output-histogram-file
  if (argc == 1)
  {
    cout << "Insufficient arguments, provide argument tree files" << endl;
    return 0;
  }
  else if (argc == 2)
  {
    afile = argv[1];
    calfile = "LabCalFileSegment.cal";
    outfile = "out.dat";
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = "out.dat";
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
  }
  else if (argc > 4)
  {
    printf("Too many arguments\n");
    return 0;
  }

  printf("Input file:%s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

  theApp=new TApplication("App", &argc, argv);

  TParserLibrary::Get()->Load();

  seg_check(afile, calfile, outfile);

  return 0;
}
