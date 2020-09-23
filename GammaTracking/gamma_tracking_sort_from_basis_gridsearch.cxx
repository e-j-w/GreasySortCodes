#include "common.h" //define all global variables here!

//Function which sorts hit positions using a pre-generated waveform basis
void sort_from_basis(const char *infile, const char *basisfileCoarse, const char *basisfileFine, const char *calfile, const char *outfile) {

  TRandom3 *rand = new TRandom3();

  TList *list = new TList;

  //read in histograms from basis file
  TFile *inp = new TFile(basisfileCoarse,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open basis file!" << endl;
    exit(-1);
  }else{
    cout << "Opened waveform basis file: " << basisfileCoarse << endl;
  }
  
  //setup histograms for the basis 
  TH1I *basisHP;
  if((basisHP = (TH1I*)inp->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the coarse waveform basis." << endl;
    exit(-1);
  }
  Int_t coarseBasisBinsR = BASIS_BINS_R_COARSE;
  Int_t coarseBasisBinsAngle = BASIS_BINS_ANGLE_COARSE;
  Int_t coarseBasisBinsZ = BASIS_BINS_Z_COARSE;
  TH1D *coarseBasis[coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ];
  for(int k = 0; k < coarseBasisBinsR; k++){
    for(int j = 0; j < coarseBasisBinsAngle; j++){
      for(int i = 0; i < coarseBasisBinsZ; i++){
        Int_t basisInd = k*coarseBasisBinsAngle*coarseBasisBinsZ + j*coarseBasisBinsZ + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/coarseBasisBinsR,(k+1)*MAX_VAL_R/coarseBasisBinsR,j*360/coarseBasisBinsAngle,(j+1)*360/coarseBasisBinsAngle,i*MAX_VAL_Z/coarseBasisBinsZ,(i+1)*MAX_VAL_Z/coarseBasisBinsZ);
        if((coarseBasis[basisInd] = (TH1D*)inp->Get(hname))==NULL){
          cout << "No coarse waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
          //cout << hname << endl;
          //getc(stdin);
        }
        //list->Add(coarseBasis[basisInd]);
      }
    }
  }
  cout << "Coarse waveform basis data read in." << endl;
  

  //read in histograms from basis file
  TFile *inp2 = new TFile(basisfileFine,"read");
  if (!inp2->IsOpen()) {
    cout << "ERROR: Could not open basis file!" << endl;
    exit(-1);
  }else{
    cout << "Opened waveform basis file: " << basisfileFine << endl;
  }
  
  //setup histograms for the basis 
  if((basisHP = (TH1I*)inp2->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the fine waveform basis." << endl;
    exit(-1);
  }
  Int_t fineBasisBinsR = BASIS_BINS_R_COARSE*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsAngle = BASIS_BINS_ANGLE_COARSE*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsZ = BASIS_BINS_Z_COARSE*FINE_BASIS_BINFACTOR;
  TH1D *fineBasis[fineBasisBinsR*fineBasisBinsAngle*fineBasisBinsZ];
  for(int k = 0; k < fineBasisBinsR; k++){
    for(int j = 0; j < fineBasisBinsAngle; j++){
      for(int i = 0; i < fineBasisBinsZ; i++){
        Int_t basisInd = k*fineBasisBinsAngle*fineBasisBinsZ + j*fineBasisBinsZ + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/fineBasisBinsR,(k+1)*MAX_VAL_R/fineBasisBinsR,j*360/fineBasisBinsAngle,(j+1)*360/fineBasisBinsAngle,i*MAX_VAL_Z/fineBasisBinsZ,(i+1)*MAX_VAL_Z/fineBasisBinsZ);
        if((fineBasis[basisInd] = (TH1D*)inp2->Get(hname))==NULL){
          cout << "No fine waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
        }
        //list->Add(fineBasis[basisInd]);
      }
    }
  }
  cout << "Fine waveform basis data read in." << endl;

  TH3D *pos3DMap = new TH3D("pos3DMap","pos3DMap",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);


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
  TTigressHit * tigress_hit;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    exit(-1);
  }

  const std::vector<Short_t> *wf, *segwf;
  bool found1, found2;
  Int_t waveform_t0;
  Double_t seg_waveform_baseline;
  Int_t numSegHits; //counter for the number of segments with a hit (ie. over the threshold energy)
  Int_t segNum;
  Double_t core_E;
  Int_t hitInd;
  Double_t wfrmSampleVal, basisSampleVal, chisq, minChisq;
  Int_t minInd, rBinCoarse, angleBinCoarse, zBinCoarse, rBinFine, angleBinFine, zBinFine;
  double rVal, angleVal, zVal;
  Int_t evtSegHP = 0; //event segment hitpattern, which will be compared against hitpatterns in the basis
  Int_t one = 1;
  //for (int jentry = 0; jentry < 100000; jentry++) {
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;
      core_E = tigress_hit->GetCharge();
      if((core_E <= 0)||(core_E > BASIS_MAX_ENERGY)) continue; //bad energy
      tigress_hit->SetWavefit();
      wf = tigress_hit->GetWaveform();
      if(wf->size()!=SAMPLES){
        cout << "Entry " << jentry << ", improper core waveform size (" << wf->size() << ")." << endl;
        continue;
      }
      TPulseAnalyzer pulse;
      pulse.SetData(*wf,0);  // Allows you to use the full TPulseAnalyzer class
      waveform_t0 = (Int_t)pulse.fit_newT0(); //in samples
      if((waveform_t0 <= 0)||(waveform_t0 >= SAMPLES-WAVEFORM_SAMPLING_WINDOW -1)){
        //this entry has an unusable risetime
        continue;
      }
      bool goodWaveforms = true;
      bool isHit = false;
      //cout << "Number of segments: " << tigress_hit->GetSegmentMultiplicity() << endl;
      numSegHits=0;
      if(tigress_hit->GetSegmentMultiplicity() == NSEG){
        //all segments have waveforms
        evtSegHP = 0;
        //check that the waveforms are the same size
        for(int i = 0; i < tigress_hit->GetSegmentMultiplicity(); i++){
          if(tigress_hit->GetSegmentHit(i).GetWaveform()->size()!=SAMPLES){
            cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
            goodWaveforms = false;
            break;
          }
          if(tigress_hit->GetSegmentHit(i).GetCharge() > SEGMENT_ENERGY_THRESHOLD){
            evtSegHP|=(one<<(tigress_hit->GetSegmentHit(i).GetSegment())); //GetSegment() is 1-indexed, evtSegHp=0 means no hits
            numSegHits++;
          }
        }
        //cout << "good: " << goodWaveforms << ", num seg hits: " << numSegHits << endl;
        if(goodWaveforms){

          minChisq = BIG_NUMBER;
          minInd = -1;

          //take action depending on the number of hit segments
          if(numSegHits==1){
            //signal best represented by a basis waveform from a position within the segment of interest
            for(int i=0; i<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i++){
              if(coarseBasis[i]!=NULL){
                //check that the event hitpattern matches the hitpattern in the basis
                //if(basisHP->GetBinContent(i+1) == evtSegHP){
                  chisq = 0;
                  for(int j=0; j<NSEG; j++){
                    segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                    segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                    seg_waveform_baseline = 0.;
                    for(int k = 0; k < BASELINE_SAMPLES; k++){
                      seg_waveform_baseline += segwf->at(k);
                    }
                    seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                    for(int k=0; k<SAMPLES; k++){
                      wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/core_E;
                      basisSampleVal = coarseBasis[i]->GetBinContent(segNum*SAMPLES + k + 1);
                      chisq += pow(fabs(wfrmSampleVal - basisSampleVal),2);
                      /*if(evtSegHP==2){
                        rBin = (Int_t)(i/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                        angleBin = (Int_t)((i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                        zBin = (i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBin*coarseBasisBinsZ);
                        cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                        cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                        cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                      }*/
                    }
                  }
                  //cout << "index " << i << ", HP: " << evtSegHP << ", chisq: " << chisq << endl;
                  if(chisq<minChisq){
                    minChisq = chisq;
                    minInd = i;
                  }
                //}
              }
            }
          }else{
            //cout << "Entry " << jentry << " contains an unhandled number of hit segments (" << numSegHits << ")." << endl;
            continue;
          }

          if(minInd >= 0){

            //a best-fit coarse basis waveform was found, figure out the corresponding bin
            rBinCoarse = (Int_t)(minInd/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
            angleBinCoarse = (Int_t)((minInd - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
            zBinCoarse = (minInd - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBinCoarse*coarseBasisBinsZ);
            //rVal = (rBin+0.5)*MAX_VAL_R/(1.0*coarseBasisBinsR);
            //angleVal = (angleBin+0.5)*360./(1.0*coarseBasisBinsAngle);
            //zVal = (zBin+0.5)*MAX_VAL_Z/(1.0*coarseBasisBinsZ);

            minInd = -1;

            if(numSegHits==1){
              //signal best represented by a basis waveform from a position within the segment of interest
              for(int i=0; i<(FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR); i++){
                rBinFine = (rBinCoarse*FINE_BASIS_BINFACTOR) + (Int_t)(i/(FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR*1.0));
                angleBinFine = (angleBinCoarse*FINE_BASIS_BINFACTOR) + (((Int_t)(i/(FINE_BASIS_BINFACTOR*1.0)) % FINE_BASIS_BINFACTOR) );
                zBinFine = (zBinCoarse*FINE_BASIS_BINFACTOR) + (i % FINE_BASIS_BINFACTOR);
                Int_t basisInd = rBinFine*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine*fineBasisBinsZ + zBinFine;
                //cout << "coarse basis bins [" << rBinCoarse << " " << angleBinCoarse << " " << zBinCoarse << "]" << endl;
                //cout << "fine basis bins [" << rBinFine << " " << angleBinFine << " " << zBinFine << "], fine basis index: " << basisInd << endl;
                if(fineBasis[basisInd]!=NULL){
                  //check that the event hitpattern matches the hitpattern in the basis
                  //if(basisHP->GetBinContent(i+1) == evtSegHP){
                    chisq = 0;
                    for(int j=0; j<NSEG; j++){
                      segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                      segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                      seg_waveform_baseline = 0.;
                      for(int k = 0; k < BASELINE_SAMPLES; k++){
                        seg_waveform_baseline += segwf->at(k);
                      }
                      seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                      for(int k=0; k<SAMPLES; k++){
                        wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/core_E;
                        basisSampleVal = fineBasis[basisInd]->GetBinContent(segNum*SAMPLES + k + 1);
                        chisq += pow(fabs(wfrmSampleVal - basisSampleVal),2);
                        /*if(evtSegHP==2){
                          rBin = (Int_t)(i/(fineBasisBinsAngle*fineBasisBinsZ*1.0));
                          angleBin = (Int_t)((i - rBin*(fineBasisBinsAngle*fineBasisBinsZ*1.0))/(fineBasisBinsZ*1.0));
                          zBin = (i - rBin*(fineBasisBinsAngle*fineBasisBinsZ) - angleBin*fineBasisBinsZ);
                          cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                          cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                          cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                        }*/
                      }
                    }
                    //cout << "index " << basisInd << ", HP: " << evtSegHP << ", chisq: " << chisq << endl;
                    if(chisq<minChisq){
                      minChisq = chisq;
                      minInd = basisInd;
                    }
                  //}
                }
              }
            }else{
              //cout << "Entry " << jentry << " contains an unhandled number of hit segments (" << numSegHits << ")." << endl;
              continue;
            }

            if(minInd >= 0){

              //a best-fit fine basis waveform was found, figure out the corresponding bin
              rBinFine = (Int_t)(minInd/(fineBasisBinsAngle*fineBasisBinsZ*1.0));
              angleBinFine = (Int_t)((minInd - rBinFine*(fineBasisBinsAngle*fineBasisBinsZ*1.0))/(fineBasisBinsZ*1.0));
              zBinFine = (minInd - rBinFine*(fineBasisBinsAngle*fineBasisBinsZ) - angleBinFine*fineBasisBinsZ);
              rVal = (rBinFine+0.5)*MAX_VAL_R/(1.0*fineBasisBinsR);
              angleVal = (angleBinFine+0.5)*360./(1.0*fineBasisBinsAngle);
              zVal = (zBinFine+0.5)*MAX_VAL_Z/(1.0*fineBasisBinsZ);
              //cout << "Filling at: x=" << rVal*cos(angleVal*M_PI/180.) << ", y=" << rVal*sin(angleVal*M_PI/180.) << ", z=" << zVal << ", minInd: " << minInd << endl;
              pos3DMap->Fill(rVal*cos(angleVal*M_PI/180.),rVal*sin(angleVal*M_PI/180.),zVal);
            }

            
          }

        }
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  inp->Close();
  inp2->Close();
}

int main(int argc, char ** argv) {

  const char *afile, *basisfileCoarse, *basisfileFine, *outfile, *calfile;

  if (argc < 2) {
    cout << endl << "This sortcode sorts interaction positions using a waveform basis (generated using the GammaTrackingMakeBasis code) using the grid search method." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasisGridSearch analysis_tree_file basis_file_coarse basis_file_fine cal_file output_file" << endl << endl;
    cout << "The analysis tree (containing the calibration data used to make the basis) is a required argument.  Omitting other arguments will cause the sortcode to fall back to default values." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    basisfileCoarse = "trackingWaveformBasisCoarse.root";
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 3) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 4) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTest.root";
  }else if (argc == 5) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = argv[4];
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 6) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = argv[4];
	  outfile = argv[5];
  } else if (argc > 6) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasis analysis_tree_file basis_file_coarse basis_file_fine cal_file output_file" << endl;
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

  cout << "Input file: " << afile << endl << "Coarse Waveform basis file: " << basisfileCoarse << endl << "Fine Waveform basis file: " << basisfileFine << endl;
  cout << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile);

  return 0;
}
