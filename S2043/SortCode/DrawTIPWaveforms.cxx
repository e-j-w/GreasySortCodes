//g++ DrawTIPWaveforms.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o DrawTIPWaveforms

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCutG.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TTip.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TEmma.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TCanvas.h"
#include "TApplication.h"

#define DrawTIP_cxx

using namespace std;

TApplication *theApp;

void SortData(char const *afile, char const *calfile, const Int_t goodPID, const Int_t channel)
{

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen())
  {
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long analentries = AnalysisTree->GetEntries();

  TTip *tip = 0;
  if (AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
  }

  //Defining Pointers
  TTipHit *tip_hit;
  const std::vector<Short_t> *wf; //for CsI waveform

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  printf("\nSorting analysis events...\n");
  for (int jentry = 0; jentry < analentries; jentry++){

    AnalysisTree->GetEntry(jentry);
    
    if(tip){
      for (int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){

        tip_hit = tip->GetTipHit(tipHitInd);
        if((channel >= 0)&&(tip_hit->GetTipChannel() != channel)){
          continue; //skip this hit
        }

        wf = tip_hit->GetWaveform();
        TPulseAnalyzer pulse;
        pulse.SetData(*wf, 0);
        if((!goodPID) || (pulse.CsIPID() > -1000.0)){
        //if((pulse.CsIPID() > 40.0) && (pulse.CsIPID() < 80.0) && (tip_hit->GetEnergy() > 2000.0)){
          if(wf->size() > 50){
            double baseline = 0;
            double max, amplitude;
            for (unsigned int i = 1; i < wf->size()-2; i++){
              if (i < 51)
                baseline += wf->at(i);
              if (i == 1)
                max = wf->at(i);
              if (max < wf->at(i))
                max = wf->at(i);
              //cout << "sample " << i << ": " << wf->at(i) << endl;
            }
            baseline = baseline / 50.0;
            amplitude = (abs(max - baseline));
            cout << "Drawing waveform" << endl << "TIP channel: " << tip_hit->GetTipChannel() << endl;
            cout << "Samples: " << wf->size() << endl << "Estimated amplitude: " << amplitude << endl;
            if((pulse.CsIPID() > -1000.0)){
              cout << "PID value: " << pulse.CsIPID() << endl;
            }else{
              cout << "PID value: " << pulse.CsIPID() << " (failed fit)" << endl;
            }
            cout << "File->Quit ROOT to cycle to the next waveform." << endl << endl;

            if((pulse.CsIPID() > -1000.0)){
              //good fit, draw the fit as well as the waveform
              cout << "PID fit parameters: " << endl;
              TCanvas *c1 = new TCanvas("c1","Histogram",200,10,1200,1000);
              //pulse.DrawT0fit();
              pulse.DrawCsIFit();
              theApp->Run(kTRUE);
              cout << endl;
            }else{
              //bad fit, draw the waveform only
              TCanvas *c1 = new TCanvas("c1","Waveform",200,10,1200,1000);
              TH1I *whist = pulse.GetWaveHist();
              char histName[256];
              whist->GetXaxis()->SetTitle("Sample Number");
              whist->SetStats(0);
              whist->Draw();
              theApp->Run(kTRUE);
            }
            
            
          }
        }

      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "No more events in this run." << endl;

}
int main(int argc, char **argv)
{

  char const *afile;
  char const *calfile;
  Int_t goodPID, channel;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0)
  {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1)
  {
    cout << "Arguments: DrawTIPWaveforms analysis_tree calibration_file good_PID channel" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    cout << "good_PID: 0 = draw all waveforms, 1 = draw only waveforms where the PID fit did not fail" << endl;
    cout << "If the channel is not specified, waveforms on all channels will be drawn." << endl;
    return 0;
  }
  else if (argc == 2)
  {
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    goodPID = 0;
    channel = -1;
    printf("Analysis file: %s\nCalibration file: %s\nTIP channel: any\n", afile, calfile);
    
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    goodPID = 0;
    channel = -1;
    printf("Analysis file: %s\nCalibration file: %s\nTIP channel: any\n", afile, calfile);
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    goodPID = atoi(argv[3]);
    channel = -1;
    printf("Analysis file: %s\nCalibration file: %s\nTIP channel: %i\n", afile, calfile, channel);
  }
  else if (argc == 5)
  {
    afile = argv[1];
    calfile = argv[2];
    goodPID = atoi(argv[3]);
    channel = atoi(argv[4]);
    printf("Analysis file: %s\nCalibration file: %s\nTIP channel: %i\n", afile, calfile, channel);
  }
  else
  {
    printf("Arguments: DrawTIPWaveforms analysis_tree calibration_file good_PID channel\n");
    return 0;
  }

  if(goodPID!=0){
    goodPID=1;
    cout << "Will only draw waveforms where PID fit doesn't fail." << endl;
  }else{
    cout << "Will draw waveforms even if PID fit fails." << endl;
  }

  theApp=new TApplication("App", &argc, argv);

  SortData(afile, calfile, goodPID, channel);

  return 0;
}
