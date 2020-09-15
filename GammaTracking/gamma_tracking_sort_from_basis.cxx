#include "common.h" //define all global variables here!

//Function which sorts hit positions using a pre-generated waveform basis
void sort_from_basis(const char *infile, const char *basisfile, const char *calfile, const char *outfile) {

  TList * list = new TList;

  //read in histograms from basis file
  TFile *inp = new TFile(basisfile,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open basis file!" << endl;
    exit(-1);
  }else{
    cout << "Opened waveform basis file: " << basisfile << endl;
  }
  
  //setup histograms for the basis 
  TH1I *basisHP;
  if((basisHP = (TH1I*)inp->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the waveform basis." << endl;
    exit(-1);
  }
  TH1D *basis[BASIS_BINS_R*BASIS_BINS_ANGLE*BASIS_BINS_Z];
  Int_t numEvtsBasis[BASIS_BINS_R*BASIS_BINS_ANGLE*BASIS_BINS_Z];
  memset(numEvtsBasis,0,sizeof(numEvtsBasis));
  for(int k = 0; k < BASIS_BINS_R; k++){
    for(int j = 0; j < BASIS_BINS_ANGLE; j++){
      for(int i = 0; i < BASIS_BINS_Z; i++){
        Int_t basisInd = k*BASIS_BINS_ANGLE*BASIS_BINS_Z + j*BASIS_BINS_Z + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/BASIS_BINS_R,(k+1)*MAX_VAL_R/BASIS_BINS_R,j*360/BASIS_BINS_ANGLE,(j+1)*360/BASIS_BINS_ANGLE,i*MAX_VAL_Z/BASIS_BINS_Z,(i+1)*MAX_VAL_Z/BASIS_BINS_Z);
        if((basis[basisInd] = (TH1D*)inp->Get(hname))==NULL){
          cout << "No waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
        }
        //list->Add(basis[basisInd]);
      }
    }
  }
  cout << "Waveform basis data read in." << endl;

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
  Double_t core_E;
  Int_t hitInd;
  Double_t wfrmSampleVal, chisq, minChisq;
  Int_t minInd, rBin, angleBin, zBin;
  double rVal, angleVal, zVal;
  Int_t evtSegHP = 0; //event segment hitpattern, which will be compared against hitpatterns in the basis
  Int_t one = 1;
  for (int jentry = 0; jentry < 100000; jentry++) {
  //for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;
      core_E = tigress_hit->GetEnergy();
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
            evtSegHP|=(one<<tigress_hit->GetSegmentHit(i).GetSegment());
            numSegHits++;
          }
        }
        if(goodWaveforms){
          
          chisq = 0.;
          minChisq = BIG_NUMBER;
          minInd = -1;

          //take action depending on the number of hit segments
          switch(numSegHits){
            case 2:
              //signal best represented by either:
              //a linear combination of 2 basis waveforms each containing a hit on one of the segments of interest
              //(for neighbouring segments) a basis waveform near the boundary of one of the two segments
              break;
            case 1:
              //signal best represented by a basis waveform from a position within the segment of interest
              for(int i=0; i<(BASIS_BINS_R*BASIS_BINS_ANGLE*BASIS_BINS_Z); i++){
                if(basis[i]!=NULL){
                  //check that the event hitpattern matches the hitpattern in the basis
                  if(basisHP->GetBinContent(i) == evtSegHP){
                    chisq=0;

                    for(int j=0; j<NSEG; j++){
                      seg_waveform_baseline = 0.;
                      for(int k = 0; k < BASELINE_SAMPLES; k++){
                        seg_waveform_baseline += tigress_hit->GetSegmentHit(j).GetWaveform()->at(k);
                      }
                      seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;

                      for(int k=BASIS_START_SAMPLE; k<SAMPLES; k++){
                        //will need to scale waveforms to match basis
                        wfrmSampleVal = (tigress_hit->GetSegmentHit(j).GetWaveform()->at(k) - seg_waveform_baseline)/core_E;
                        chisq += pow(wfrmSampleVal - basis[i]->GetBinContent((tigress_hit->GetSegmentHit(j).GetSegment())*SAMPLES + k),2);
                      }
                    }
                    //cout << "index " << i << ", chisq: " << chisq << endl;
                    if(chisq<minChisq){
                      minChisq = chisq;
                      minInd = i;
                    }
                  }
                }
                  
              }
              break;
            case 0:
              //do nothing
              //cout << "Entry " << jentry << " contains no hit segments." << endl;
              break;
            default:
              cout << "Entry " << jentry << " contains an unhandled number of hit segments (" << numSegHits << ")." << endl;
              break;
          }

          if(minInd >= 0){

            //a best-fit basis waveform was found, figure out the corresponding hit location
            rBin = minInd/(BASIS_BINS_ANGLE*BASIS_BINS_Z);
            angleBin = (minInd - rBin*(BASIS_BINS_ANGLE*BASIS_BINS_Z))/BASIS_BINS_Z;
            zBin = (minInd - rBin*(BASIS_BINS_ANGLE*BASIS_BINS_Z) - angleBin*BASIS_BINS_Z);
            rVal = (rBin+0.5)*MAX_VAL_R/(1.0*BASIS_BINS_R);
            angleVal = (angleBin+0.5)*360./(1.0*BASIS_BINS_ANGLE);
            zVal = (zBin+0.5)*MAX_VAL_Z/(1.0*BASIS_BINS_Z);
            //cout << "Filling at: x=" << rVal*cos(angleVal*M_PI/180.) << ", y=" << rVal*sin(angleVal*M_PI/180.) << ", z=" << zVal << ", minInd: " << minInd << endl;
            pos3DMap->Fill(rVal*cos(angleVal*M_PI/180.),rVal*sin(angleVal*M_PI/180.),zVal);
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
}

int main(int argc, char ** argv) {

  const char *afile, *basisfile, *outfile, *calfile;

  if (argc < 2) {
    cout << endl << "This sortcode sorts interaction positions using a waveform basis (generated using the GammaTrackingMakeBasis code)." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasis analysis_tree_file basis_file cal_file output_file" << endl << endl;
    cout << "The analysis tree (containing the calibration data used to make the basis) is a required argument.  Omitting other arguments will cause the sortcode to fall back to default values." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    basisfile = "trackingWaveformBasis.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 3) {
	  afile = argv[1];
    basisfile = argv[2];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 4) {
	  afile = argv[1];
    basisfile = argv[2];
	  calfile = argv[3];
	  outfile = "trackingBasisSortTest.root";
  } else if (argc == 5) {
	  afile = argv[1];
    basisfile = argv[2];
	  calfile = argv[3];
	  outfile = argv[4];
  } else if (argc > 5) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasis analysis_tree_file basis_file cal_file output_file" << endl;
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

  cout << "Input file: " << afile << endl << "Waveform basis file: " << basisfile << endl << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  sort_from_basis(afile, basisfile, calfile, outfile);

  return 0;
}
