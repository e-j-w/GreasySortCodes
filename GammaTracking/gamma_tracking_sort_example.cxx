#include "common.h" //define all global variables here!

//function which generates a mapping between ordering parameters and real spatial coordinates
//and saves this mapping to disk
void sort_test(const char *infile, const char *mapfile, const char *calfile, const char *outfile) {

  TList * list = new TList;

  TH3D *pos3DMap = new TH3D("pos3DMap","pos3DMap",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);

  //read in histograms from map file
  TFile *inp = new TFile(mapfile,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open map file!" << endl;
    exit(-1);
  }
  char hname[20];
  TH1 *rMap[NSEG], *angleMap[NSEG*MAX_VAL_R/BIN_WIDTH_R], *zMap[NSEG*(MAX_VAL_R/BIN_WIDTH_R)*(MAX_VAL_ANGLE/BIN_WIDTH_ANGLE)];
  for(int k = 0; k < NSEG; k++){
    sprintf(hname,"rMapSeg%i",k);
    if((rMap[k] = (TH1*)inp->Get(hname))==NULL){
      cout << "No r coordinate map for segment " << k << endl;
    }
    for(int j = 0; j < MAX_VAL_R/BIN_WIDTH_R; j++){
      sprintf(hname,"angleMapSeg%ir%ito%i",k,j*BIN_WIDTH_R,(j+1)*BIN_WIDTH_R);
      if((angleMap[k*MAX_VAL_R/BIN_WIDTH_R + j] = (TH1*)inp->Get(hname))==NULL){
        cout << "No angle coordinate map for segment " << k << ", radial bin " << j << endl;
      }
      for(int i = 0; i < MAX_VAL_ANGLE/BIN_WIDTH_ANGLE; i++){
        sprintf(hname,"zMapSeg%ir%ito%iangle%ito%i",k,j*BIN_WIDTH_R,(j+1)*BIN_WIDTH_R,i*BIN_WIDTH_ANGLE,(i+1)*BIN_WIDTH_ANGLE);
        if((zMap[k*(MAX_VAL_ANGLE/BIN_WIDTH_ANGLE)*(MAX_VAL_R/BIN_WIDTH_R) + j*MAX_VAL_ANGLE/BIN_WIDTH_ANGLE + i] = (TH1*)inp->Get(hname))==NULL){
          cout << "No z coordinate map for segment " << k << ", radial bin " << j << ", angle bin " << k << endl;
        }
      }
    } 
  }

  cout << "Map file data read in." << endl;

  //setup histograms for the mapped parameters
  TH1D *rMappedHist[NSEG], *angleMappedHist[NSEG], *zMappedHist[NSEG];
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
  TTigressHit * tigress_hit;
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

  const std::vector<Short_t> *wf;
  bool found1, found2;
  Int_t waveform_t0;
  Int_t one;
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

        Int_t posNum = tigress_hit->GetDetector()-1;
        Int_t coreNum = tigress_hit->GetCrystal();
        Int_t segNum = tigress_hit->GetSegmentHit(i).GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

        //calculate all ordering parameters (see ordering_parameter_calc.cxx)
        double rho = calc_ordering(tigress_hit,i,jentry,waveform_t0,0);
        if(rho == BAD_RETURN){
          continue;
        }
        double phi = calc_ordering(tigress_hit,i,jentry,waveform_t0,1);
        if(phi == BAD_RETURN){
          continue;
        }
        double zeta = calc_ordering(tigress_hit,i,jentry,waveform_t0,2);
        if(zeta == BAD_RETURN){
          continue;
        }

        map_hit_counter++;

        //here is where the mapping happens
        //cout << "rho: " << rho << ", phi: " << phi << ", zeta: " << zeta << endl;
        //cout << "map bins: " << rMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(rho) << ", " << angleMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(phi) << ", " << zMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(zeta) << endl;
        double r=-1.;
        double angle=-1.;
        double z=-1.;
        if(rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          r = rMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(rMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(rho));
        }
        if(angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          angle = angleMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(angleMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(phi));
        }
        if(zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]!=NULL){
          z = zMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->GetBinContent(zMap[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->FindBin(zeta));
        }
        if((r>=0.)&&(z>=0.)&&(angle>=0.)){
          rMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(r);
          angleMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(angle);
          zMappedHist[NCORE*NSEG*posNum + NSEG*coreNum + segNum]->Fill(z);
          angle += 90.0*(segNum%4);
          //if(z<10)
            pos3DMap->Fill(r*cos(angle*M_PI/180.),r*sin(angle*M_PI/180.),z);
          //cout << "r: " << r << ", angle: " << angle << ", z: " << z << endl;
          /*if(r==0.){
            cout << "rho: " << rho << endl;
          }
          cout << "x: " << r*cos(angle*M_PI/180.) << ", y: " << r*sin(angle*M_PI/180.) << ", z: " << z << endl;*/
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
