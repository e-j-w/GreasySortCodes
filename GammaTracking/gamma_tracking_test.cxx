//g++ wf_sing.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o WFSing

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TGriffin.h"
#include "TGRSIDetectorHit.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TPulseAnalyzer.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

#define     NPOS   16 //number of positions in the array
#define     NCORE  4  //number of cores per position
#define     NSEG   8  //number of segments per core

//lists of adjacent segments in the TIGRESS array (zero-indexed)
Int_t phiAdjSeg1[8] = {3,0,1,2,7,4,5,6};
Int_t phiAdjSeg2[8] = {1,2,3,0,5,6,7,4};
Int_t    zAdjSeg[8] = {4,5,6,7,0,1,2,3};

//function which generates a mapping between ordering parameters and real spatial coordinates
//and saves this mapping to disk
void generate_mapping(char const * infile, char const * calfile, char const * outfile) {

  TList * list = new TList;

  //setup density functions for the real spatial coordinates r, angle, and z in cylindrical coordinates
  //these are rough guesses, and would be better replaced by simulated distributions
  //for r, assume proportional to r^2 (crystal radius 30mm)
  TH1D *rDistHist = new TH1D("radius distribution","radius distribution",6,0,30);
  for(int i=0;i<6;i++){
    double rVal = rDistHist->GetBinCenter(i+1);
    rDistHist->SetBinContent(i+1,rVal*rVal);
  }
  rDistHist->Scale(1.0/rDistHist->Integral());
  list->Add(rDistHist);
  //for angle, assume flat distribution
  TH1D *angleDistHist = new TH1D("angle distribution","angle distribution",6,0,90);
  for(int i=0;i<6;i++){
    angleDistHist->SetBinContent(i+1,1.0);
  }
  angleDistHist->Scale(1.0/angleDistHist->Integral());
  list->Add(angleDistHist);
  //for z, assume exponential, decreasing with depth (what is the decay constant?)
  //lateral segmentation is 30mm from crystal front, crystals are 90mm long, 
  //so need 2 distributions corresponding to front and back segments
  TH1D *zDistHistFront = new TH1D("z distribution front","z distribution front",6,0,30);
  for(int i=0;i<6;i++){
    double decConst = 0.01;
    double zVal = zDistHistFront->GetBinCenter(i+1);
    zDistHistFront->SetBinContent(i+1,exp(-decConst*zVal));
  }
  zDistHistFront->Scale(1.0/zDistHistFront->Integral());
  list->Add(zDistHistFront);
  TH1D *zDistHistBack = new TH1D("z distribution back","z distribution back",6,30,90);
  for(int i=0;i<6;i++){
    double decConst = 0.01;
    double zVal = zDistHistBack->GetBinCenter(i+1);
    zDistHistBack->SetBinContent(i+1,exp(-decConst*zVal));
  }
  zDistHistBack->Scale(1.0/zDistHistBack->Integral());
  list->Add(zDistHistBack);
  TH1 *rDistHistC = rDistHist->GetCumulative();
  TH1 *angleDistHistC = angleDistHist->GetCumulative();
  TH1 *zDistHistFrontC = zDistHistFront->GetCumulative();
  TH1 *zDistHistBackC = zDistHistBack->GetCumulative();
  list->Add(rDistHistC);
  list->Add(angleDistHistC);
  list->Add(zDistHistFrontC);
  list->Add(zDistHistBackC);

  //setup histograms for the ordering parameters
  //ROOT histograms are used to store this data since ROOT provides useful
  //methods such as GetCumulative() and GetQuantiles() which will be used later
  TH1D *rhoHist[NPOS*NCORE*NSEG], *phiHist[NPOS*NCORE*NSEG], *zetaHist[NPOS*NCORE*NSEG];
  char hname[20];
  for(int i = 0; i < NPOS; i++){
    for(int j = 0; j < NCORE; j++){
      for(int k = 0; k < NSEG; k++){
        sprintf(hname,"rhoPos%iCore%iSeg%i",i,j,k);
        rhoHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("rhoPos%iCore%iSeg%i",i,j,k),2048,-2048,2048);
        sprintf(hname,"phiPos%iCore%iSeg%i",i,j,k);
        phiHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("phiPos%iCore%iSeg%i",i,j,k),2048,-0.1,0.1);
        sprintf(hname,"zetaPos%iCore%iSeg%i",i,j,k);
        zetaHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("zetaPos%iCore%iSeg%i",i,j,k),2048,-0.1,0.1);
        /*list->Add(rhoHist[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(phiHist[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(zetaHist[NCORE*NSEG*i + NSEG*j + k]);*/
      }
    } 
  }

  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    printf("Opening file failed, aborting\n");
    return;
  }
  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  printf("%i tree files, details:\n", AnalysisTree->GetNtrees());
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TTigress * tigress = 0;
  TTigressHit * tigress_hit, * tigress_hit2;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    return;
  }

  Int_t samples = 100; //number of samples per waveform
  Int_t sampling_window = 10; //number of waveform samples used to construct ordering parameters

  const std::vector<Short_t> *wf, *segwf, *segwf2, *segwf3;
  bool found1, found2;
  Int_t waveform_t0;
  Int_t one;
  Int_t offset = 0;
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (one = 0; one < tigress->GetMultiplicity(); one++) {
      tigress_hit = tigress->GetTigressHit(one);
      if(tigress_hit->GetKValue() != 700) continue;
      tigress_hit->SetWavefit();
      wf = tigress_hit->GetWaveform();
      TPulseAnalyzer pulse;
      pulse.SetData(*wf,0);  // Allows you to use the full TPulseAnalyzer class
      waveform_t0 = (Int_t)pulse.fit_newT0(); //in samples
      if((waveform_t0 <= 0)||(waveform_t0 >= samples-sampling_window-1)){
        //this entry has an unusable risetime
        continue;
      }
      for(int i = 0; i < tigress_hit->GetSegmentMultiplicity(); i++)
      {
        TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(i);
        Int_t posNum = tigress_hit->GetArrayNumber();
        Int_t coreNum = tigress_hit->GetCrystal();
        Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed
        //cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
        segwf = segment_hit.GetWaveform();

        //construct rho, the ordering parameter for the radius
        //see Eq. 4 of NIM A 729 (2013) 198-206
        double sampleAvg = 0.;
        double rho = 0.;
        double dno = 0.; //placeholder for denominator
        for(int j=0;j<sampling_window;j++){
          sampleAvg += waveform_t0+j;
          dno += (segwf->at(waveform_t0+j+1) - segwf->at(waveform_t0+j-1))/2.0;
        }
        sampleAvg /= sampling_window*1.0;
        for(int j=0;j<sampling_window;j++){
          rho += pow(waveform_t0+j - sampleAvg,3.0)*(segwf->at(waveform_t0+j+1) - segwf->at(waveform_t0+j-1))/2.0;
        }
        rho /= dno;

        //contruct phi, the ordering parameter for the angle
        //see Eq. 3 of NIM A 729 (2013) 198-206
        double phi = 0.;
        found1 = false;
        found2 = false;
        for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
          if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
            found1=true;
            segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
          }
          if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
            found2=true;
            segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }
        }
        if((!found1)||(!found2)){
          cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute phi parameter." << endl;
          cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
          continue;
        }
        for(int j=0;j<sampling_window;j++){
          phi += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
          dno += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
        }
        phi /= dno;

        //contruct zeta, the ordering parameter for the z direction
        //see Eq. 2 of NIM A 729 (2013) 198-206 (modified here)
        double zeta = 0.;
        found1 = false;
        for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
          if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
            found1=true;
            segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }
        }
        if(!found1){
          cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
          cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
          continue;
        }
        for(int j=0;j<sampling_window;j++){
          zeta += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) - segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
          dno += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) + segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
        }
        zeta /= dno;

        rhoHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(rho);
        phiHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(phi);
        zetaHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(zeta);

      }
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  //generate normalized cumulative distributions of all ordering parameters
  cout << "Generating cumulative distributions of ordering parameters..." << endl;
  TH1 *rhoHistC[NPOS*NCORE*NSEG], *phiHistC[NPOS*NCORE*NSEG], *zetaHistC[NPOS*NCORE*NSEG];
  for(int i = 0; i < NPOS; i++){
    for(int j = 0; j < NCORE; j++){
      for(int k = 0; k < NSEG; k++){
        rhoHist[NCORE*NSEG*i + NSEG*j + k]->Scale(1.0/rhoHist[NCORE*NSEG*i + NSEG*j + k]->Integral());
        rhoHistC[NCORE*NSEG*i + NSEG*j + k] = rhoHist[NCORE*NSEG*i + NSEG*j + k]->GetCumulative();
        phiHist[NCORE*NSEG*i + NSEG*j + k]->Scale(1.0/phiHist[NCORE*NSEG*i + NSEG*j + k]->Integral());
        phiHistC[NCORE*NSEG*i + NSEG*j + k] = phiHist[NCORE*NSEG*i + NSEG*j + k]->GetCumulative();
        zetaHist[NCORE*NSEG*i + NSEG*j + k]->Scale(1.0/zetaHist[NCORE*NSEG*i + NSEG*j + k]->Integral());
        zetaHistC[NCORE*NSEG*i + NSEG*j + k] = zetaHist[NCORE*NSEG*i + NSEG*j + k]->GetCumulative();
        list->Add(rhoHistC[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(phiHistC[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(zetaHistC[NCORE*NSEG*i + NSEG*j + k]);
      }
    } 
  }

  //generate ordering parameter to spatial coordinate maps
  cout << "Generating ordering parameter to spatial coordinate maps..." << endl;
  Double_t *qVal;
  TH1 *rMap[NPOS*NCORE*NSEG], *angleMap[NPOS*NCORE*NSEG], *zMap[NPOS*NCORE*NSEG];
  for(int i = 0; i < NPOS; i++){
    for(int j = 0; j < NCORE; j++){
      for(int k = 0; k < NSEG; k++){
        sprintf(hname,"rMapPos%iCore%iSeg%i",i,j,k);
        rMap[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("rMapPos%iCore%iSeg%i",i,j,k),6,0,30);
        sprintf(hname,"angleMapPos%iCore%iSeg%i",i,j,k);
        angleMap[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("angleMapPos%iCore%iSeg%i",i,j,k),6,0,90);
        sprintf(hname,"zMapPos%iCore%iSeg%i",i,j,k);
        zMap[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("zMapPos%iCore%iSeg%i",i,j,k),6,0,90);
        list->Add(rMap[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(angleMap[NCORE*NSEG*i + NSEG*j + k]);
        list->Add(zMap[NCORE*NSEG*i + NSEG*j + k]);
      }
    } 
  }

  cout << "Histograms written, sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();
}

int main(int argc, char ** argv) {

  char const * afile;
  char const * outfile;
  char const * calfile;
  cout << "Starting sortcode..." << endl;

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0) {
	  grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
	  cout << "Insufficient arguments, provide analysis tree files" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
	  calfile = "LabCalFileSegment.cal";
	  outfile = "trackingTest.root";
  } else if (argc == 3) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = "trackingTest.root";
  } else if (argc == 4) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = argv[3];
  } else if (argc > 4) {
	  printf("Too many arguments\n");
	  return 0;
  }

  printf("Input file:%s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

  TParserLibrary::Get()->Load();

  generate_mapping(afile, calfile, outfile);

  return 0;
}
