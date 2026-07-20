//g++ seg_check_drawwaveforms.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -o SegCheckDraw
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

void seg_check(char const *infile, char const *calfile, char const *outfile, int wfsample)
{

  TList *list = new TList;

  // Energy Spectra
  TH1F *gE = new TH1F("gE", "gE", 2048, 0, 2048);
  list->Add(gE);
  TH1D *gSeg[10];
  char shname[20];
  for (int iii = 0; iii < 8; iii++)
  {
    sprintf(shname, "Segment%d", iii);
    gSeg[iii] = new TH1D(shname, Form("Segment d%1.1i", iii), 2048, 0, 2048);
    list->Add(gSeg[iii]);
  }

  // Waveform Spectra
  TH2F *enVamp = new TH2F("enVamp", " ", 2048, 0, 2048, 2048, 0, 2048);
  list->Add(enVamp);
  TH2F *segVamp[8];
  char hname[20];
  for (int iii = 0; iii < 8; iii++)
  {
    sprintf(hname, "segVa%i", iii);
    segVamp[iii] = new TH2F(hname, Form("Segment %1.1i", iii + 1), 2048, 0, 2048, 2048, 0, 2048);
    list->Add(segVamp[iii]);
    segVamp[iii]->GetXaxis()->SetTitle("WaveForm Fit");
    segVamp[iii]->GetYaxis()->SetTitle("Evaluated Energy");
  }

  TList *wavelist = new TList;
  TH1F *wave[100];
  for (int iii = 0; iii < 100; iii++)
  {
    sprintf(hname, "h%d", iii);
    wave[iii] = new TH1F(hname, Form("WaveForm %1.1i", iii), wfsample, 0, wfsample);
    wavelist->Add(wave[iii]);
  }

  TH1F *segwave[8][100];
  for (int jjj = 0; jjj < 8; jjj++)
  {
    for (int iii = 0; iii < 100; iii++)
    {
      sprintf(hname, "h%d_%d", jjj, iii);
      segwave[jjj][iii] = new TH1F(hname, Form("Segment WaveForm %1.1i %1.1i", jjj + 1, iii), wfsample, 0, wfsample);
      wavelist->Add(segwave[jjj][iii]);
    }
  }

  TH1F *segmult = new TH1F("segmult", "", 30, 0, 30);
  list->Add(segmult);
  TH2F *segsamp = new TH2F("segsamp", "", 8, 0, 8, 750, 0, 750);
  list->Add(segsamp);

  TH2F *tdiffT0 = new TH2F("tdiffT0", "", 10, 0, 10, 2000, -1000, 1000);
  list->Add(tdiffT0);
  TH2F *tdiffCfd = new TH2F("tdiffCfd", "", 10, 0, 10, 2000, -1000, 1000);
  list->Add(tdiffCfd);
  TH2F *tdiffCfdT0 = new TH2F("tdiffCfdt0", "", 10, 0, 10, 2000, -1000, 1000);
  list->Add(tdiffCfdT0);
  TH2F *tdiff = new TH2F("tdiff", "", 10, 0, 10, 2000, -1000, 1000);
  list->Add(tdiff);
  TH2F *tdiffTS = new TH2F("tdiffTS", "", 10, 0, 10, 2000, -1000, 1000);
  list->Add(tdiffTS);

  TFile *inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen())
  {
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
  if (AnalysisTree->FindBranch("TTigress"))
  {
    AnalysisTree->SetBranchAddress("TTigress", &tigr);
  }
  else
  {
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
  for (int jentry = 0; jentry < (nentries - 1); jentry++)
  {
    tree->GetEntry(jentry);
    for (int one = 0; one < tigr->GetMultiplicity(); one++)
    {
      TCanvas *c1 = new TCanvas("c1","Histogram",200,10,1200,1000);
      c1->Divide(3,3);

      tigress_hit = tigr->GetTigressHit(one);
      gE->Fill(tigress_hit->GetEnergy());
      tigress_hit->SetWavefit();
      wf = tigress_hit->GetWaveform();
      baseline = 0;
      for (unsigned int i = 0; i < wf->size(); i++)
      {
        if (i < 20)
          baseline += wf->at(i);
        if (i == 0)
          max = wf->at(i);
        if (max < wf->at(i))
          max = wf->at(i);
        if (waveform_counter < 100)
          wave[waveform_counter]->Fill(i, wf->at(i));
      }
      if(wf->size() <= 0){ continue; }
      amplitude = (abs(max - baseline * 0.05));
      printf("\nCore\n----\n");
      printf("FV waveform size: %i\nAmplitude: %f",wf->size(),amplitude);
      /*for (unsigned int i = 0; i < wf->size(); i++){
        printf("sample %i: %i\n",i,wf->at(i));
      }*/
      enVamp->Fill(amplitude, tigress_hit->GetEnergy());
      segmult->Fill(tigress_hit->GetSegmentMultiplicity());
      TPulseAnalyzer pulse;
      pulse.SetData(*wf, 0);
      //tmpTime = pulse.fit_newT0() * 10.;
      //if(tmpTime <= 0.) continue;
      tmpCfd = tigress_hit->GetTime();
      tmpCfdT0 = tigress_hit->GetTime() - ((((tigress_hit->GetTimeStamp()) + gRandom->Uniform()) * tigress_hit->GetTimeStampUnit()) - 10. * 50.); //offset by pretrigger

      c1->cd(1);
      pulse.DrawT0fit();
      
      //printf("Fit t0: %f ns\nCFD t0: %f ns\nCFD time: %f ns\nTimestamp time: %f ns\n",tmpTime,tmpCfdT0,tigress_hit->GetTime(),(tigress_hit->GetTimeStamp()+0.0)*tigress_hit->GetTimeStampUnit());
      printf("\nSegments\n--------\n");
      /*TH1I *cwhist = pulse.GetWaveHist();
      cwhist->SetStats(0);
      cwhist->SetTitle("Core");
      cwhist->GetXaxis()->SetTitle("Sample Number");
      cwhist->Draw();*/
      if(tigress_hit->GetSegmentMultiplicity()!=8) continue;

      for (int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++)
      {
        const TDetectorHit* segment_hit = &tigress_hit->GetSegmentHit(j);
        gSeg[segment_hit->GetSegment() - 1]->Fill(segment_hit->GetEnergy());
        segwf = segment_hit->GetWaveform();
        baseline = 0;
        segsamp->Fill(segment_hit->GetSegment() - 1, segwf->size());
        for (unsigned int i = 0; i < segwf->size(); i++)
        {
          if (i < 20)
            baseline += segwf->at(i);
          if (i == 0)
            max = segwf->at(i);
          if (max < segwf->at(i))
            max = segwf->at(i);
          if (segwaveform_counter[segment_hit->GetSegment() - 1] < 100)
            segwave[segment_hit->GetSegment() - 1][segwaveform_counter[segment_hit->GetSegment() - 1]]->Fill(i, segwf->at(i));
        }
        amplitude = (abs(max - baseline * 0.05));
        printf("segment waveform size: %i, amplitude: %f\n",segwf->size(),amplitude);
        /*for (unsigned int i = 0; i < segwf->size(); i++){
          printf("sample %i: %i\n",i,segwf->at(i));
        }*/
        //	if(segwaveform_counter[segment_hit->GetSegment()-1] < 100) cout << segwaveform_counter[segment_hit->GetSegment()-1] << "\t" << segment_hit->GetSegment()-1 << "\t" << amplitude << "\t" << segment_hit->GetEnergy() << endl;
        segVamp[segment_hit->GetSegment() - 1]->Fill(amplitude, segment_hit->GetEnergy());
        TPulseAnalyzer segpulse;
        segpulse.SetData(*segwf, 0);
        tmpTimeSeg = segpulse.fit_newT0() * 10.;
        tmpCfdSeg = segment_hit->GetTime();
        tmpCfdSegT0 = segment_hit->GetTime() - ((((segment_hit->GetTimeStamp()) + gRandom->Uniform()) * segment_hit->GetTimeStampUnit()) - 10. * 50.); //offset by pretrigger
        c1->cd(2+j);
        TH1I *whist = segpulse.GetWaveHist();
        if(whist != NULL){
          snprintf(histName,255,"Segment %i",segment_hit->GetSegment());
          whist->SetTitle(histName);
          whist->GetXaxis()->SetTitle("Sample Number");
          whist->SetStats(0);
          whist->Draw();
        }
        //printf("Segment %i fit t0: %f ns\nSegment %i CFD t0: %f ns\nSegment %i CFD time: %f ns\nSegment %i timestamp time: %f ns\n",segment_hit->GetSegment(),tmpTimeSeg,segment_hit->GetSegment(),tmpCfdSegT0,segment_hit->GetSegment(),segment_hit->GetTime(),segment_hit->GetSegment(),(segment_hit->GetTimeStamp()+0.0)*segment_hit->GetTimeStampUnit());

        if (tigress_hit->GetEnergy() > 100 && segment_hit->GetEnergy() > 100)
        //if (tigress_hit->GetEnergy() > 100 && amplitude > 400)
        {
          tdiffCfd->Fill(segment_hit->GetSegment(), tmpCfd - tmpCfdSeg);
          tdiffCfdT0->Fill(segment_hit->GetSegment(), tmpCfdT0 - tmpCfdSegT0);
          tdiffT0->Fill(segment_hit->GetSegment(), tmpTime - tmpTimeSeg);
          tdiff->Fill(segment_hit->GetSegment(), (tmpTime - tmpCfd) - (tmpTimeSeg - tmpCfdSeg));
          tdiffTS->Fill(segment_hit->GetSegment(), tigress_hit->GetTimeStampNs() - segment_hit->GetTimeStampNs());
        }
        segwaveform_counter[segment_hit->GetSegment() - 1]++;
      }
      waveform_counter++;

      if(tigress_hit->GetEnergy() > 100)
        if(tigress_hit->GetSegmentMultiplicity()>0)
          theApp->Run(kTRUE);
      //getc(stdin);

    }
    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "TIGRESS Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete"
           << "\r" << flush;
  }
  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";
  cout << "Histograms written, sorting complete" << endl;

  cout << gE->Integral(100, -1) << endl;
  for (int i = 0; i < 8; i++)
    cout << gSeg[i]->Integral(100, -1) << endl;
  cout << "Writing histograms to " << outfile << endl;
  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  TDirectory *wavedir = myfile->mkdir("Waveforms");
  wavedir->cd();
  wavelist->Write();
  myfile->cd();
  myfile->Close();
}

int main(int argc, char **argv)
{

  char const *afile;
  char const *outfile;
  char const *calfile;
  int samples;
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
    outfile = "segments.root";
    samples = 100;
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = "segments.root";
    samples = 100;
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    samples = 100;
  }
  else if (argc == 5)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    samples = atoi(argv[4]);
  }
  else if (argc > 5)
  {
    printf("Too many arguments\n");
    return 0;
  }

  printf("Input file:%s\nCalibration file: %s\nOutput file: %s\nNum Samples: %d\n", afile, calfile, outfile, samples);

  theApp=new TApplication("App", &argc, argv);

  TParserLibrary::Get()->Load();

  seg_check(afile, calfile, outfile, samples);

  return 0;
}
