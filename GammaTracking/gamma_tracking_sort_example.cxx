#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TLeaf.h"
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

#define     N_BINS_ORDERING 2048 //number of bins to use when discretizing ordering parameter

//lists of adjacent segments in the TIGRESS array (zero-indexed)
Int_t phiAdjSeg1[8] = {3,0,1,2,7,4,5,6};
Int_t phiAdjSeg2[8] = {1,2,3,0,5,6,7,4};
Int_t    zAdjSeg[8] = {4,5,6,7,0,1,2,3};

//function which generates a mapping between ordering parameters and real spatial coordinates
//and saves this mapping to disk
void sort_test(const char *infile, const char *mapfile, const char *calfile, const char *outfile) {

  TList * list = new TList;

  //read in histograms from map file
  TFile *inp = new TFile(mapfile,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open map file!" << endl;
    exit(-1);
  }
  char hname[20];
  TH1 *rMap[NPOS*NCORE*NSEG], *angleMap[NPOS*NCORE*NSEG], *zMap[NPOS*NCORE*NSEG];
  for(int i = 0; i < NPOS; i++){
    for(int j = 0; j < NCORE; j++){
      for(int k = 0; k < NSEG; k++){
        sprintf(hname,"rMapPos%iCore%iSeg%i",i,j,k);
        if((rMap[NCORE*NSEG*i + NSEG*j + k] = (TH1*)inp->Get(hname))==NULL){
          cout << "No r coordinate map for position " << i << ", core " << j << ", seg " << k << endl;
        }
        sprintf(hname,"angleMapPos%iCore%iSeg%i",i,j,k);
        if((angleMap[NCORE*NSEG*i + NSEG*j + k] = (TH1*)inp->Get(hname))==NULL){
          cout << "No angle coordinate map for position " << i << ", core " << j << ", seg " << k << endl;
        }
        sprintf(hname,"zMapPos%iCore%iSeg%i",i,j,k);
        if((zMap[NCORE*NSEG*i + NSEG*j + k] = (TH1*)inp->Get(hname))==NULL){
          cout << "No z coordinate map for position " << i << ", core " << j << ", seg " << k << endl;
        }
      }
    } 
  }

  cout << "Map file data read in." << endl;

  //setup histograms for the mapped parameters
  TH1D *rMappedHist[NPOS*NCORE*NSEG], *angleMappedHist[NPOS*NCORE*NSEG], *zMappedHist[NPOS*NCORE*NSEG];
  for(int i = 0; i < NPOS; i++){
    for(int j = 0; j < NCORE; j++){
      for(int k = 0; k < NSEG; k++){
        if(rMap[NCORE*NSEG*i + NSEG*j + k]!=NULL){
          sprintf(hname,"rMappedPos%iCore%iSeg%i",i,j,k);
          rMappedHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("rMappedPos%iCore%iSeg%i",i,j,k),40,0,40);
          list->Add(rMappedHist[NCORE*NSEG*i + NSEG*j + k]);
        }
        if(angleMap[NCORE*NSEG*i + NSEG*j + k]!=NULL){
          sprintf(hname,"angleMappedPos%iCore%iSeg%i",i,j,k);
          angleMappedHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("angleMappedPos%iCore%iSeg%i",i,j,k),30,0,90);
          list->Add(angleMappedHist[NCORE*NSEG*i + NSEG*j + k]);
        }
        if(zMap[NCORE*NSEG*i + NSEG*j + k]!=NULL){
          sprintf(hname,"zMappedPos%iCore%iSeg%i",i,j,k);
          zMappedHist[NCORE*NSEG*i + NSEG*j + k] = new TH1D(hname,Form("zMappedPos%iCore%iSeg%i",i,j,k),90,0,90);
          list->Add(zMappedHist[NCORE*NSEG*i + NSEG*j + k]);
        }
      }
    } 
  }

  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    cout << "ERROR: Could not open analysis tree file!" << endl;
    exit(-1);
  }
  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  cout << AnalysisTree->GetNtrees() << " tree files, details:" << endl;
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  cout << "Reading calibration file: " << calfile << endl;
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TTigress * tigress = 0;
  TTigressHit * tigress_hit, * tigress_hit2;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    exit(-1);
  }

  Int_t samples = 100; //number of samples per waveform
  Int_t sampling_window = 10; //number of waveform samples used to construct ordering parameters

  Int_t hit_counter = 0;
  Int_t map_hit_counter = 0;

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
      samples = wf->size();
      TPulseAnalyzer pulse;
      pulse.SetData(*wf,0);  // Allows you to use the full TPulseAnalyzer class
      waveform_t0 = (Int_t)pulse.fit_newT0(); //in samples
      if((waveform_t0 <= 0)||(waveform_t0 >= samples-sampling_window-1)){
        //this entry has an unusable risetime
        continue;
      }
      for(int i = 0; i < tigress_hit->GetSegmentMultiplicity(); i++)
      {

        hit_counter++;
        TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(i);

        Int_t posNum = tigress_hit->GetDetector()-1;
        Int_t coreNum = tigress_hit->GetCrystal();
        Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

        //cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
        segwf = segment_hit.GetWaveform();
        if((posNum < 0)||(posNum > 15)){
          cout << "Entry " << jentry << ", invalid array position: " << posNum << endl;
          continue;
        }
        if(segwf->size() != samples){
          cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
          continue;
        }

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
        if(rho!=rho){
          cout << "Entry " << jentry << ", cannot compute rho parameter (NaN)." << endl;
          continue;
        }

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
        }else if((segwf->size() != samples)||(segwf2->size() != samples)){
          cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
          continue;
        }
        for(int j=0;j<sampling_window;j++){
          phi += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
          dno += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
        }
        phi /= dno;
        if(phi!=phi){
          cout << "Entry " << jentry << ", cannot compute phi parameter (NaN)." << endl;
          continue;
        }

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
        }else if(segwf3->size() != samples){
          cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
          continue;
        }
        for(int j=0;j<sampling_window;j++){
          zeta += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) - segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
          dno += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) + segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
        }
        zeta /= dno;
        if(zeta!=zeta){
          cout << "Entry " << jentry << ", cannot compute zeta parameter (NaN)." << endl;
          continue;
        }

        if(segNum>3){
          //back segment, reverse sign to make zeta increase with z
          zeta *= -1.;
        }

        map_hit_counter++;

        //here is where the mapping happens
        if(rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindFirstBinAbove(rho)));
        }
        if(angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindFirstBinAbove(phi)));
        }
        if(zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindFirstBinAbove(zeta)));
        }
      }
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << map_hit_counter << " of " << hit_counter << " hits retained (" << 100*map_hit_counter/hit_counter << " %)." << endl;

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  inp->Close();
}

int main(int argc, char ** argv) {

  const char *afile, *mapfile, *outfile, *calfile;

  // Input-chain-file, output-histogram-file
  if (argc < 2) {
    cout << endl << "This sortcode sorts resconstructed hit positions, using the map file generated using the GammaTrackingMakeMap code." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortExample analysis_tree_file map_file cal_file output_file" << endl;
    cout << "The analysis tree is a required argument.  Omitting other arguments will cause the sortcode to fall back to default values." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    mapfile = "trackingMap.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingSortTest.root";
  } else if (argc == 3) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingSortTest.root";
  } else if (argc == 4) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = argv[3];
	  outfile = "trackingSortTest.root";
  } else if (argc == 5) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = argv[3];
	  outfile = argv[4];
  } else if (argc > 5) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingSortExample analysis_tree_file map_file cal_file output_file" << endl;
	  return 0;
  }

  cout << "Starting sortcode..." << endl;

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0) {
	  grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  cout << "Input file: " << afile << endl << "Simulation data file: " << mapfile << endl << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  sort_test(afile, mapfile, calfile, outfile);

  return 0;
}
