//#include "TCanvas.h"
//#include "TApplication.h"
#include "GammaTrackingTIGRESS.h" //define all global variables here!

//TApplication *theApp;

TH3D *pos3DMap, *pos3DMapClover, *pos3DMapAbs;
TH2D *posXYMapBottom, *posXYMapTop, *posXZMap, *posYZMap;
TH2D *posXYMapBottom137Cs, *posXYMapTop137Cs, *posXZMap137Cs, *posYZMap137Cs;
TH1D *tigE, *dopplerESeg, *dopplerEGT;
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

  TVector3 recoil_vec;
  recoil_vec.SetX(0.); recoil_vec.SetY(0.); recoil_vec.SetZ(1.);
  double beta = 0.05;

  //for (int jentry = 0; jentry < 100000; jentry++) {
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;
      hit_counter++;
      TVector3 posVec = GT_get_pos_gridsearch(tigress_hit,&trackingBasis);
      if(posVec.X()!=BAD_RETURN){
        //cout << "Filling: " << posVec.X() << " " << posVec.Y() << " " << posVec.Z() << endl;
        pos3DMap->Fill(posVec.X(),posVec.Y(),posVec.Z());
        TVector3 posVecAbs = GT_transform_position_to_absolute(tigress_hit,&posVec);
        pos3DMapAbs->Fill(posVecAbs.X(),posVecAbs.Y(),posVecAbs.Z());
        TVector3 posVecClover = GT_transform_position_to_clover(tigress_hit,&posVec);
        pos3DMapClover->Fill(posVecClover.X(),posVecClover.Y(),posVecClover.Z());
        if(posVec.Z() <= 30.){
          posXYMapBottom->Fill(posVec.X(),posVec.Y());
        }else{
          posXYMapTop->Fill(posVec.X(),posVec.Y());
        }
        posXZMap->Fill(posVec.X(),posVec.Z());
        posYZMap->Fill(posVec.Y(),posVec.Z());
        //map positions of hits corresponding to 137Cs photopeak
        if((tigress_hit->GetEnergy() > 659)&&(tigress_hit->GetEnergy() < 664)){
          if(posVec.Z() <= 30.){
            posXYMapBottom137Cs->Fill(posVec.X(),posVec.Y());
          }else{
            posXYMapTop137Cs->Fill(posVec.X(),posVec.Y());
          }
          posXZMap137Cs->Fill(posVec.X(),posVec.Z());
          posYZMap137Cs->Fill(posVec.Y(),posVec.Z());
        }
        double egt = GT_get_doppler(beta,&recoil_vec,tigress_hit,&posVec);
        if(egt>5.){
          dopplerEGT->Fill(egt);
        }
        sort_hit_counter++;
      }
      double eseg = tigress_hit->GetDoppler(beta,&recoil_vec);
      if(eseg>5.){
        dopplerESeg->Fill(eseg);
        tigE->Fill(tigress_hit->GetEnergy());
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << sort_hit_counter << " of " << hit_counter << " hits retained (" << 100*sort_hit_counter/hit_counter << " %)." << endl;
}

//Function which sorts hit positions using a pre-generated waveform basis
void sort_from_basis(const char *infile, const char *basisfileCoarse, const char *basisfileFine, const char *calfile, const char *outfile, bool inpList) {

  //read in histograms from basis files
  TFile *coarseBasisInp = new TFile(basisfileCoarse,"read");
  TFile *fineBasisInp = new TFile(basisfileFine,"read");
  GT_import_basis(coarseBasisInp,fineBasisInp,&trackingBasis);

  //setup output histograms
  TList *list = new TList;
  pos3DMap = new TH3D("pos3DMap","Core Position Map",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);
  pos3DMapClover = new TH3D("pos3DMapClover","Clover Positon Map",80,-80,80,80,-80,80,40,-10,100);
  list->Add(pos3DMapClover);
  pos3DMapAbs = new TH3D("pos3DMapAbsolute","Absolute Positon Map",40,-300,300,40,-300,300,40,-300,300);
  list->Add(pos3DMapAbs);
  posXYMapBottom = new TH2D("posXYMapSeg0-3","posXYMapSeg0-3",40,-40,40,40,-40,40);
  list->Add(posXYMapBottom);
  posXYMapTop = new TH2D("posXYMapSeg4-7","posXYMapSeg4-7",40,-40,40,40,-40,40);
  list->Add(posXYMapTop);
  posXZMap = new TH2D("posXZMap","posXZMap",40,-40,90,40,-10,100);
  list->Add(posXZMap);
  posYZMap = new TH2D("posYZMap","posYZMap",40,-40,90,40,-10,100);
  list->Add(posYZMap);
  posXYMapBottom137Cs = new TH2D("posXYMapSeg0-3 - 137Cs full energy","posXYMapSeg0-3 - 137Cs full energy",40,-40,40,40,-40,40);
  list->Add(posXYMapBottom137Cs);
  posXYMapTop137Cs = new TH2D("posXYMapSeg4-7 - 137Cs full energy","posXYMapSeg4-7 - 137Cs full energy",40,-40,40,40,-40,40);
  list->Add(posXYMapTop137Cs);
  posXZMap137Cs = new TH2D("posXZMap - 137Cs full energy","posXZMap - 137Cs full energy",40,-40,90,40,-10,100);
  list->Add(posXZMap137Cs);
  posYZMap137Cs = new TH2D("posYZMap - 137Cs full energy","posYZMap - 137Cs full energy",40,-40,90,40,-10,100);
  list->Add(posYZMap137Cs);
  tigE = new TH1D("tigE","TIGRESS energy (uncorrected)",16384,0,4096);
  list->Add(tigE);
  dopplerESeg = new TH1D("dopplerESeg","Doppler Corrected Energy (segment position)",16384,0,4096);
  list->Add(dopplerESeg);
  dopplerEGT = new TH1D("dopplerEGT","Doppler Corrected Energy (tracked position)",16384,0,4096);
  list->Add(dopplerEGT);

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
  coarseBasisInp->Close();
  fineBasisInp->Close();

}

int main(int argc, char ** argv) {

  const char *afile, *basisfileCoarse, *basisfileFine, *outfile, *calfile;
  char *ext;

  if (argc < 2) {
    cout << endl << "This sortcode sorts interaction positions using a waveform basis (generated using the GammaTrackingMakeBasis code) using the grid search method." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasisGridSearch analysis_tree_file basis_file_coarse basis_file_fine cal_file output_file" << endl << endl;
    cout << "The analysis tree (containing the experimental data used to be sorted) is a required argument.  Omitting other arguments will cause the sortcode to fall ";
    cout << "back to default values.  The analysis tree can be a single ROOT file, or a list of ROOT files (using file extension .list) can be specified instead." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    basisfileCoarse = "trackingWaveformBasisCoarse.root";
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  } else if (argc == 3) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  } else if (argc == 4) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  }else if (argc == 5) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = argv[4];
	  outfile = "trackingBasisSortTestGridSearch.root";
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

  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/

  //theApp=new TApplication("App", &argc, argv);
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, false);
  }

  return 0;
}
