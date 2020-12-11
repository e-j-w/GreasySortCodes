//#include "TCanvas.h"
//#include "TApplication.h"
#include "GammaTrackingTIGRESS.h" //define all global variables here!

//TApplication *theApp;

TH2D *posXYMap, *posXYMapGridSearch, *posXYMapDirect;
TH3D *pos3DMap, *pos3DMapGridSearch, *pos3DMapDirect;
GT_map trackingMap;
GT_basis trackingBasis;

void sortData(TFile *inputfile, const char *calfile){

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

  Int_t hit_counter = 0;
  Int_t sort_hit_counter = 0;
  Int_t sort_hit_counter_gridsearch = 0;
  Int_t sort_hit_counter_direct = 0;

  //for (int jentry = 0; jentry < 100000; jentry++) {
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;
      hit_counter++;
      TVector3 posVec = GT_get_pos_gridsearch(tigress_hit,&trackingBasis); //sort using grid search
      if(posVec.X()!=BAD_RETURN){
        pos3DMap->Fill(posVec.X(),posVec.Y(),posVec.Z());
        posXYMap->Fill(posVec.X(),posVec.Y());
        pos3DMapGridSearch->Fill(posVec.X(),posVec.Y(),posVec.Z());
        posXYMapGridSearch->Fill(posVec.X(),posVec.Y());
        sort_hit_counter++;
        sort_hit_counter_gridsearch++;
      }else{
        posVec = GT_get_pos_direct(tigress_hit,&trackingMap); //sort using direct method
        if(posVec.X()!=BAD_RETURN){
          pos3DMap->Fill(posVec.X(),posVec.Y(),posVec.Z());
          posXYMap->Fill(posVec.X(),posVec.Y());
          pos3DMapDirect->Fill(posVec.X(),posVec.Y(),posVec.Z());
          posXYMapDirect->Fill(posVec.X(),posVec.Y());
          sort_hit_counter++;
          sort_hit_counter_direct++;
        }
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << sort_hit_counter << " of " << hit_counter << " hits retained (" << 100*sort_hit_counter/hit_counter << " %)." << endl;
  cout << sort_hit_counter_gridsearch << " hits sorted using grid search (" << 100*sort_hit_counter_gridsearch/sort_hit_counter << " % of sorted hits)." << endl;
  cout << sort_hit_counter_direct << " hits sorted using direct method (" << 100*sort_hit_counter_direct/sort_hit_counter << " % of sorted hits)." << endl;
}

//Function which sorts hit positions using a pre-generated waveform basis
void sort_hybrid(const char *infile, const char *mapfile, const char *basisfileCoarse, const char *basisfileFine, const char *calfile, const char *outfile, bool inpList) {

  //read in histograms from map file
  TFile *mapInp = new TFile(mapfile,"read");
  GT_import_map(mapInp,&trackingMap);

  //read in histograms from basis files
  TFile *coarseBasisInp = new TFile(basisfileCoarse,"read");
  TFile *fineBasisInp = new TFile(basisfileFine,"read");
  GT_import_basis(coarseBasisInp,fineBasisInp,&trackingBasis);

  //setup output histograms
  TList *list = new TList;
  pos3DMap = new TH3D("pos3DMap","pos3DMap",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);
  pos3DMapGridSearch = new TH3D("pos3DMapGridSearch","pos3DMapGridSearch",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMapGridSearch);
  pos3DMapDirect = new TH3D("pos3DMapDirect","pos3DMapDirect",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMapDirect);
  posXYMap = new TH2D("posXYMap","posXYMap",40,-40,40,40,-40,40);
  list->Add(posXYMap);
  posXYMapGridSearch = new TH2D("posXYMapGridSearch","posXYMapGridSearch",40,-40,40,40,-40,40);
  list->Add(posXYMapGridSearch);
  posXYMapDirect = new TH2D("posXYMapDirect","posXYMapDirect",40,-40,40,40,-40,40);
  list->Add(posXYMapDirect);

  //sort data from individual analysis tree or list of trees
  if(inpList){
    FILE *listFile;
    char name[256];
    if((listFile=fopen(infile,"r"))==NULL){
      cout << "ERROR: Could not open analysis tree list file!" << endl;
      exit(-1);
    }
    while(fscanf(listFile,"%s",name)!=EOF){
      TFile * inputfile = new TFile(name, "READ");
      if (!inputfile->IsOpen()) {
        cout << "ERROR: Could not open analysis tree file (" << name << ")" << endl;
        exit(-1);
      }
      sortData(inputfile, calfile);
      inputfile->Close();
    }
    fclose(listFile);
  }else{
    TFile * inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      cout << "ERROR: Could not open analysis tree file (" << infile << ")" << endl;
      exit(-1);
    }
    sortData(inputfile, calfile);
    inputfile->Close();
  }

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  mapInp->Close();
  coarseBasisInp->Close();
  fineBasisInp->Close();

}

int main(int argc, char ** argv) {

  const char *afile, *mapfile, *basisfileCoarse, *basisfileFine, *outfile, *calfile;
  char *ext;

  if (argc < 2) {
    cout << endl << "This sortcode sorts interaction positions using a waveform basis (generated using the GammaTrackingMakeBasis code) using the grid search method, ";
    cout << "and falls back to the direct method if the drif search fails." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortHybrid analysis_tree_file map_file basis_file_coarse basis_file_fine cal_file output_file" << endl << endl;
    cout << "The analysis tree (containing the experimental data used to be sorted) is a required argument.  Omitting other arguments will cause the sortcode to fall ";
    cout << "back to default values.  The analysis tree can be a single ROOT file, or a list of ROOT files (using file extension .list) can be specified instead." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    mapfile = "trackingMap.root";
    basisfileCoarse = "trackingWaveformBasisCoarse.root";
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestHybrid.root";
  } else if (argc == 3) {
	  afile = argv[1];
    mapfile = argv[2];
    basisfileCoarse = "trackingWaveformBasisCoarse.root";
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestHybrid.root";
  } else if (argc == 4) {
	  afile = argv[1];
    mapfile = argv[2];
    basisfileCoarse = argv[3];
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestHybrid.root";
  }else if (argc == 5) {
	  afile = argv[1];
    mapfile = argv[2];
    basisfileCoarse = argv[3];
    basisfileFine = argv[4];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestHybrid.root";
  } else if (argc == 6) {
	  afile = argv[1];
    mapfile = argv[2];
    basisfileCoarse = argv[3];
    basisfileFine = argv[4];
	  calfile = argv[5];
	  outfile = "trackingBasisSortTestHybrid.root";
  } else if (argc == 7) {
	  afile = argv[1];
    mapfile = argv[2];
    basisfileCoarse = argv[3];
    basisfileFine = argv[4];
	  calfile = argv[5];
	  outfile = argv[6];
  } else if (argc > 7) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingSortHybrid analysis_tree_file map_file basis_file_coarse basis_file_fine cal_file output_file" << endl;
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

  cout << "Input file: " << afile << endl << "Map data file: " << mapfile << endl << "Coarse Waveform basis file: " << basisfileCoarse << endl;
  cout << "Fine Waveform basis file: " << basisfileFine << endl;
  cout << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/

  //theApp=new TApplication("App", &argc, argv);
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    sort_hybrid(afile, mapfile, basisfileCoarse, basisfileFine, calfile, outfile, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    sort_hybrid(afile, mapfile, basisfileCoarse, basisfileFine, calfile, outfile, false);
  }

  return 0;
}
