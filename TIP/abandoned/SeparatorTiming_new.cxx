//Sort code to separate data based on TIP and TIGRESS timing
//windows are defined in common.h

#define SeparatorTiming_cxx
#include "common.h"
#include "SeparatorTiming.h"

using namespace std;

//verify that a separated tree contains entries matching an original tree
void SeparatorTiming::VerifyTimingSep(TTree *sepTree, TTree *origTree){

  cout << "Verifying trees..." << endl;

  TTigress *tigress = 0;
  TTigress *tigressSep = 0;
  if(origTree->FindBranch("TTigress")){
    origTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "ERROR: Branch 'TTigress' not found in original tree, could not verify tree!" << endl;
    return;
  }
  if(sepTree->FindBranch("TTigress")){
    sepTree->SetBranchAddress("TTigress", &tigressSep);
  }else{
    cout << "ERROR: Branch 'TTigress' not found in separated tree, could not verify tree!" << endl;
    return;
  }

  TTip *tip = 0;
  TTip *tipSep = 0;
  if(origTree->FindBranch("TTip")){
    origTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "ERROR: Branch 'TTip' not found in original tree, could not verify tree!" << endl;
    return;
  }
  if(sepTree->FindBranch("TTip")){
    sepTree->SetBranchAddress("TTip", &tipSep);
  }else{
    cout << "ERROR: Branch 'TTip' not found in separated tree, could not verify tree!" << endl;
    return;
  }

  Long64_t sepEntry = 0;

  for(Long64_t jentry = 0; jentry < origTree->GetEntries(); jentry++){

    if(origTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not verify entry " << jentry << " in original tree." << endl;
      getc(stdin);
      continue;
    }else{
      if(tigress && tip){
        if(tip->GetMultiplicity()>MAXNUMTIPHIT){
          cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
          continue;
        }
        if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
          cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
          continue;
        }

        uint64_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup
        if(passedtimeGate&(1ULL<<TIPTIGFLAG)){
          if(sepTree->GetEntry(sepEntry) == 0){
            //entry not read successfully
            cout << "WARNING: could not verify entry " << sepEntry << " in separated tree." << endl;
            getc(stdin);
          }else{

            if(tigressSep && tipSep){
              //tigressSep->ClearTransients();
              //tigressSep->ResetAddback();
              //tigressSep->BuildHits();
              bool pass = false;
              if(tigress->GetAddbackMultiplicity() == tigressSep->GetAddbackMultiplicity()){
                if(tip->GetMultiplicity() == tipSep->GetMultiplicity()){
                  pass = true;
                  for(int i=0;i<tip->GetMultiplicity();i++){
                    if(tip->GetTipHit(i)->GetTimeStamp() != tipSep->GetTipHit(i)->GetTimeStamp()){
                      pass = false;
                      break;
                    }
                  }
                  for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
                    if(tigress->GetAddbackHit(i)->GetTimeStamp() != tigressSep->GetAddbackHit(i)->GetTimeStamp()){
                      pass = false;
                      break;
                    }
                  }
                }
              }

              /*if(sepEntry==31)
                pass = false;*/

              if(pass == false){
                cout << "Original event " << jentry << " and separated event " << sepEntry << " do not match!" << endl;
                for(int i = 0; i < tip->GetMultiplicity() ; i++){
                  cout << " Original TIP hit " << i << ", ts: " << tip->GetTipHit(i)->GetTimeStamp() << ", t: " << tip->GetTipHit(i)->GetTime() 
                  << ", fit t: " << getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger) << ", K: " << tip->GetTipHit(i)->GetKValue() << ", E: " << tip->GetTipHit(i)->GetEnergy() << ", ch: " << tip->GetTipHit(i)->GetTipChannel() << endl;
                }
                for(int i = 0; i < tigress->GetMultiplicity(); i++){
                  cout << " Original TIGRESS hit " << i << ", ts: " << tigress->GetTigressHit(i)->GetTimeStamp() << ", K: " 
                  << tigress->GetTigressHit(i)->GetKValue() << ", BGO: " << tigress->GetTigressHit(i)->BGOFired() << ", chan: " << tigress->GetTigressHit(i)->GetChannel() << ", E: " << tigress->GetTigressHit(i)->GetEnergy() 
                  << ", ANum: " << tigress->GetTigressHit(i)->GetArrayNumber() << ", name: " << tigress->GetTigressHit(i)->GetName() << endl;
                  for(int j = 0; j < tigress->GetTigressHit(i)->GetNSegments(); j++){
                    cout << "  Segment hit " << j << ", name: " << tigress->GetTigressHit(i)->GetSegmentHit(j).GetName() << ", E: " << tigress->GetTigressHit(i)->GetSegmentHit(j).GetEnergy() << endl; 
                  }
                }
                for(int i = 0; i < tigress->GetAddbackMultiplicity(); i++){
                  cout << " Original TIGRESS AB hit " << i << ", ts: " << tigress->GetAddbackHit(i)->GetTimeStamp() << ", t: " << tigress->GetAddbackHit(i)->GetTime() << ", K: " 
                  << tigress->GetAddbackHit(i)->GetKValue() << ", BGO: " << tigress->GetAddbackHit(i)->BGOFired() << ", E: " << tigress->GetAddbackHit(i)->GetEnergy() 
                  << ", ANum: " << tigress->GetAddbackHit(i)->GetArrayNumber() << endl;
                }
                cout << endl;
                for(int i = 0; i < tipSep->GetMultiplicity(); i++){
                  cout << " Separated TIP hit " << i << ", ts: " << tipSep->GetTipHit(i)->GetTimeStamp() << ", t: " << tipSep->GetTipHit(i)->GetTime() 
                  << ", fit t: " << getTipFitTime(tipSep->GetTipHit(i),tip_waveform_pretrigger) << ", K: " << tipSep->GetTipHit(i)->GetKValue() << ", E: " << tipSep->GetTipHit(i)->GetEnergy() << ", ch: " << tipSep->GetTipHit(i)->GetTipChannel() << endl;
                }
                for(int i = 0; i < tigressSep->GetMultiplicity(); i++){
                  cout << " Separated TIGRESS hit " << i << ", ts: " << tigressSep->GetTigressHit(i)->GetTimeStamp() << ", K: " 
                  << tigressSep->GetTigressHit(i)->GetKValue() << ", BGO: " << tigressSep->GetTigressHit(i)->BGOFired() << ", chan: " << tigressSep->GetTigressHit(i)->GetChannel() << ", E: " << tigressSep->GetTigressHit(i)->GetEnergy() 
                  << ", ANum: " << tigressSep->GetTigressHit(i)->GetArrayNumber() << ", name: " << tigressSep->GetTigressHit(i)->GetName() << endl;
                  for(int j = 0; j < tigressSep->GetTigressHit(i)->GetNSegments(); j++){
                    cout << "  Segment hit " << j << ", name: " << tigressSep->GetTigressHit(i)->GetSegmentHit(j).GetName() << ", E: " << tigressSep->GetTigressHit(i)->GetSegmentHit(j).GetEnergy() << endl; 
                  }
                  //tigressSep->GetTigressHit(i)->Print();
                }
                for(int i = 0; i < tigressSep->GetAddbackMultiplicity(); i++){
                  cout << " Separated TIGRESS AB hit " << i << ", ts: " << tigressSep->GetAddbackHit(i)->GetTimeStamp() << ", t: " << tigressSep->GetAddbackHit(i)->GetTime() << ", K: " 
                  << tigressSep->GetAddbackHit(i)->GetKValue() << ", BGO: " << tigressSep->GetAddbackHit(i)->BGOFired() << ", E: " << tigressSep->GetAddbackHit(i)->GetEnergy() 
                  << ", ANum: " << tigressSep->GetAddbackHit(i)->GetArrayNumber() << endl;
                }
                getc(stdin);
              }
            }else{
              cout << "WARNING: No TIGRESS or TIP data in entry " << sepEntry << " of the separated tree!" << endl;
            }
          }
          sepEntry++;
        }
      }else{
        cout << "WARNING: No TIGRESS or TIP data in entry " << jentry << " of the original tree!" << endl;
      }

    }
    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << origTree->GetEntries() << ", " << 100 * jentry / origTree->GetEntries() << "% complete" << "\r" << flush;
  }

  cout << "Verified " << sepEntry << " entries in the separated tree event-by-event against the original tree." << endl;

}

void SeparatorTiming::SortData(char const *afile, char const *calfile, char const *outfile){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  std::map<TClass*, TDetector**> det_map;
  TChain *inpChain = (TChain *)analysisfile->Get("AnalysisTree");
  long int totalEntries = inpChain->GetEntries();
  TObjArray *array = inpChain->GetListOfBranches();
  for(int x=0;x<array->GetSize();x++) {
    TBranch *b = (TBranch*)array->At(x);
    if(b) {
      TClass *c = TClass::GetClass(b->GetName());
      if(c) {
        cout << "Found " << b->GetName() << endl;
        TDetector** det = new TDetector*;
        *det = NULL;
        det_map[c] = det;
        inpChain->SetBranchAddress(b->GetName(),det_map[c]);
      }
    }
  }

  //setup the output file
  //TFile *myfile = new TFile(outfile, "RECREATE");
  //setup the output tree
  //TTree *outTree = AnalysisTree->CloneTree(0);
  //outTree->CopyAddresses(AnalysisTree,kTRUE);
  /*TTree *outTree = new TTree("AnalysisTree", "AnalysisTree"); //renaming tree or branches doesn't help
  TTigress *tigressSep = 0;
  TTip *tipSep = 0;
  outTree->Branch("TTigress","TTigress",&tigressSep,320000); //defining/filling only tigress doesn't help
  outTree->Branch("TTip","TTip",&tipSep,320000);*/

  //TEntryList *sepList = new TEntryList("EntryList",afile,AnalysisTree);
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  Long64_t numSeparatedEvents = 0;

  printf("\nSorting analysis events...\n");
  for(Long64_t entriesRead = 0; entriesRead < totalEntries; entriesRead++){

    if(inpChain->GetEntry(entriesRead) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << entriesRead << endl; 
      continue;
    }

    std::shared_ptr<TUnpackedEvent> event = std::make_shared<TUnpackedEvent>();
    for(auto& elem : det_map) {
      TDetector* det = *elem.second;
      /*
      //if(!det->TestBit(TDetector::kUnbuilt)){
      if(det->IsBuilt()) {
        //event->AddDetector(det);
        std::cout << "HERE" << std::endl;
        //event->AddDetector(std::shared_ptr<TDetector>(det));
      } 
      else {
        
        //if(det->Timestamp()!=-1 && det->Size()!=0) {
        //std::cout << det->IsA()->GetName() << " was not present in this event (TS="
        //	  << det->Timestamp() << ")"
        //	  << std::endl;
        //}
        
        //std::cout << "I don't know what this condition means..." << std::endl;
        delete det;
      }
      */
      event->AddDetector(std::shared_ptr<TDetector>(det));
    }

    if(tigress && tip){

      if(tip->GetMultiplicity()>MAXNUMTIPHIT){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }else if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }else{

        uint64_t passedtimeGate = passesTimeGate(tigress,tip); //also rejects pileup
        if(passedtimeGate&(1ULL<<TIPTIGFLAG)){
          
          //tigress->ClearTransients();
          //tigress->ResetAddback();
          //tigress->BuildHits();

          //cout << "addr: " << tigress->GetTigressHit(0)->GetChannel() << endl;

          /*for(int i = 0; i < tigress->GetMultiplicity(); i++){
            tigress->GetTigressHit(i)->ClearChannel();
            tigress->GetTigressHit(i)->ClearEnergy();
          }*/
          
          /*tigress->Copy(*tigressSep);
          tip->Copy(*tipSep);*/
          
          
          //cout << "orig addr: " << tigress->GetTigressHit(0)->GetChannel() << endl;
          //cout << "sep addr: " << tigressSep->GetTigressHit(0)->GetChannel() << endl;

          //tigressSep->GetTigressHit(0)->SetEnergy(tigress->GetTigressHit(0)->GetEnergy());
          //tigress->BuildHits();
          //tigressSep->BuildHits();
          //tipSep->BuildHits();
          
          //cout << "Filling..." << endl;
          outTree->Fill(); //from all indications the original and separated data are the same at this point, but the segment data in the output tree is incorrectly filled
          //the problem persists with the latest ROOT version (6.26/04) built with C++14
          //cout << "Done filling." << endl;

          /*cout << "orig hit: " << endl;
          tigress->GetTigressHit(0)->Print();
          cout << "sep hit: " << endl;
          tigressSep->GetTigressHit(0)->Print();
          getc(stdin);*/

          //sepList->Enter(jentry,AnalysisTree);
          
          numSeparatedEvents++;
          
        }
      }

    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  cout << "Number of separated events: " << numSeparatedEvents << endl;
  cout << "Entries written in separated output data: " << outTree->GetEntries() << endl;
  if(numSeparatedEvents != outTree->GetEntries()){
    cout << "\nNumber of written entires not equal!\n\n" << endl;
    getc(stdin);
  }

  cout << "Writing histograms to " << outfile << endl;

  myfile->Write("",TObject::kOverwrite);
  myfile->Close();
  analysisfile->Close();

  //verify the number of events actually written
  TFile *checkfile = new TFile(outfile, "READ"); //Opens Analysis Tree
  if (!checkfile->IsOpen()){
    printf("Reopening file %s failed, aborting\n", afile);
    return;
  }
  TTree *sepTree = (TTree *)checkfile->Get("AnalysisTree");
  cout << "Verified number of written entries: " << sepTree->GetEntries() << endl;
  if(numSeparatedEvents != sepTree->GetEntries()){
    cout << "\nNumber of written entires not equal!\n\n" << endl;
    getc(stdin);
  }

  TFile *origfile = new TFile(afile, "READ"); //Opens Analysis Tree
  if (!origfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }
  TTree *unsepTree = (TTree *)origfile->Get("AnalysisTree");

  VerifyTimingSep(sepTree,unsepTree);

  checkfile->Close();
  origfile->Close();

}
int main(int argc, char **argv){

  SeparatorTiming *mysort = new SeparatorTiming();

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Arguments: SeparatorTiming analysis_tree calibration_file output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else if (argc == 2){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    
  }else if (argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if (argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else{
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
