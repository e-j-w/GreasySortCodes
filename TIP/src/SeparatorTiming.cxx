//Sort code to separate data based on TIP and TIGRESS timing
//windows are defined in common.h

#define SeparatorTiming_cxx
#include "common.h"
#include "evt_fmt.h"
#include "SeparatorTiming.h"

using namespace std;


uint64_t *mapping;

//verify that a separated tree contains entries matching an original tree
void SeparatorTiming::VerifyTimingSep(TTree *origTree, const char *sepfile, const char *calfile){

  cout << "Verifying separated data..." << endl;

  sorted_evt sortedEvt;
  uint8_t footerVal;
  FILE *inp = fopen(sepfile, "rb");

  TTigress *tigress = 0;
  if(origTree->FindBranch("TTigress")){
    origTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "ERROR: Branch 'TTigress' not found in original tree, could not verify tree!" << endl;
    return;
  }

  TTip *tip = 0;
  if(origTree->FindBranch("TTip")){
    origTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "ERROR: Branch 'TTip' not found in original tree, could not verify tree!" << endl;
    return;
  }

  TTigressHit *add_hit;
  TTipHit *tip_hit;

  //read the expected number of events in the separated data
  uint64_t numSepEvts = 0U;
  fread(&numSepEvts,sizeof(uint64_t),1,inp);
  uint64_t sepEntry = 0;

  for(uint64_t jentry = 0; jentry < numSepEvts; jentry++){

    if(origTree->GetEntry(mapping[jentry]) == 0){
      //entry not read successfully
      cout << "WARNING: could not verify entry " << jentry << " in original tree." << endl;
      continue;
    }

    memset(&sortedEvt,0,sizeof(sorted_evt));
    footerVal = 0;
    fread(&sortedEvt.header,sizeof(evt_header),1,inp);
    for(int i = 0; i<sortedEvt.header.numTigABHits;i++){
      fread(&sortedEvt.tigHit[i],sizeof(tigab_hit),1,inp);
    }
    for(int i = 0; i<sortedEvt.header.numCsIHits;i++){
      fread(&sortedEvt.csiHit[i],sizeof(csi_hit),1,inp);
    }
    fread(&footerVal,sizeof(uint8_t),1,inp);
    if(footerVal != 227U){
      printf("ERROR: invalid footer value in separated event %lu (%u)!\n", sepEntry, footerVal);
      cout << "Corresponding original entry: " << mapping[jentry] << endl;
      printf("Original addback multiplicity: %u\n", (uint8_t)tigress->GetAddbackMultiplicity());
      printf("Original TIP multiplicity: %u\n", (uint8_t)tip->GetMultiplicity());
      printf("Separated addback multiplicity: %u\n", sortedEvt.header.numTigABHits);
      printf("Separated TIP multiplicity: %u\n", sortedEvt.header.numCsIHits);
      exit(-1);
    }
    /*if(((uint8_t)tigress->GetAddbackMultiplicity() != sortedEvt.header.numTigABHits)||((uint8_t)tip->GetMultiplicity() != sortedEvt.header.numCsIHits)){
      printf("ERROR: invalid hit multiplicity in separated event %lu!\n", sepEntry);
      cout << "Corresponding original entry: " << mapping[jentry] << endl;
      printf("Original addback multiplicity: %u\n", (uint8_t)tigress->GetAddbackMultiplicity());
      printf("Original TIP multiplicity: %u\n", (uint8_t)tip->GetMultiplicity());
      printf("Separated addback multiplicity: %u\n", sortedEvt.header.numTigABHits);
      printf("Separated TIP multiplicity: %u\n", sortedEvt.header.numCsIHits);
      exit(-1);
    }*/

    //check event
    for(int i = 0; i<tigress->GetAddbackMultiplicity();i++){
      add_hit = tigress->GetAddbackHit(i);
      if(floorf(add_hit->GetEnergy()) != floorf(sortedEvt.tigHit[i].energy)){
        cout << "Entry " << mapping[jentry] << " addback hit " << i << " energy mismatch: [" << add_hit->GetEnergy() << "," << sortedEvt.tigHit[i].energy << "]" << endl;
      }
      if(fabs(add_hit->GetTime() - sortedEvt.tigHit[i].timeNs) > 1.0){
        cout << "Entry " << mapping[jentry] << " addback hit " << i << " time mismatch: [" << add_hit->GetTime() << "," << sortedEvt.tigHit[i].timeNs << "]" << endl;
      }
    }

    for(int i = 0; i<tip->GetMultiplicity();i++){
      tip_hit = tip->GetTipHit(i);
      if(floorf(tip_hit->GetEnergy()) != floorf(sortedEvt.csiHit[i].energy)){
        cout << "Entry " << mapping[jentry] << " CsI hit " << i << " energy mismatch: [" << tip_hit->GetEnergy() << "," << sortedEvt.csiHit[i].energy << "]" << endl;
      }
      if(floor(getTipFitTime(tip_hit,tip_waveform_pretrigger)) != floor(sortedEvt.csiHit[i].timeNs)){
        cout << "Entry " << mapping[jentry] << " CsI hit " << i << " time mismatch: [" << getTipFitTime(tip_hit,tip_waveform_pretrigger) << "," << sortedEvt.csiHit[i].timeNs << "]" << endl;
      }
    }

    sepEntry++;
      
    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Verified entry " << mapping[jentry] << " of " << origTree->GetEntries() << ", " << 100 * mapping[jentry] / origTree->GetEntries() << "% complete" << "\r" << flush;
  }

  fclose(inp);

  if(sepEntry != numSepEvts){
    cout << "ERROR: number of separated events (" << sepEntry << ") is inconsistent with expected value (" << numSepEvts << ")." << endl;
    exit(-1);
  }

  cout << "Verified " << sepEntry << " entries in the separated data event-by-event against the original tree." << endl;

}

void SeparatorTiming::SortData(const char *afile, const char *calfile, const char *outfile){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TTree *AnalysisTree = (TTree *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();
  mapping = (uint64_t*)malloc(analentries*sizeof(uint64_t));

  TTigress *tigress = 0;
  if(AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "ERROR: Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
    return;
  }

  TTip *tip = 0;
  if(AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "ERROR: Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
    return;
  }

  TTigressHit *add_hit;
  TTipHit *tip_hit;

  //setup the output file
  FILE *out = fopen(outfile, "wb");

  sorted_evt sortedEvt;
  uint8_t footerVal = 227U;
  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  TChannel::ReadCalFile(calfile);

  uint64_t numSeparatedEvents = 0;

  printf("\nSorting analysis events...\n");
  for(uint64_t jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    if(tigress && tip){

      if((tip->GetMultiplicity()>MAXNUMTIPHIT)||(tip->GetMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }else if((tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT)||(tigress->GetAddbackMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }else{

        uint64_t passedtimeGate = passesTimeGate(tigress,tip,1,2); //also rejects pileup
        if(passedtimeGate&(1ULL<<TIPTIGFLAG)){

          mapping[numSeparatedEvents] = jentry;

          memset(&sortedEvt,0,sizeof(sorted_evt));
          sortedEvt.header.numTigABHits = (uint8_t)tigress->GetAddbackMultiplicity();
          sortedEvt.header.numCsIHits = (uint8_t)tip->GetMultiplicity();
          sortedEvt.header.numRFHits = (uint8_t)0;

          fwrite(&sortedEvt.header,sizeof(evt_header),1,out);
          
          for(int i = 0; i<sortedEvt.header.numTigABHits;i++){
            add_hit = tigress->GetAddbackHit(i);
            sortedEvt.tigHit[i].energy = (float)add_hit->GetEnergy();
            sortedEvt.tigHit[i].timeNs = (double)add_hit->GetTime();
            sortedEvt.tigHit[i].core = (uint8_t)add_hit->GetArrayNumber();
            sortedEvt.tigHit[i].seg = (uint8_t)add_hit->GetFirstSeg();
            fwrite(&sortedEvt.tigHit[i],sizeof(tigab_hit),1,out);
          }

          for(int i = 0; i<sortedEvt.header.numCsIHits;i++){
            tip_hit = tip->GetTipHit(i);
            sortedEvt.csiHit[i].energy = (float)tip_hit->GetEnergy();
            sortedEvt.csiHit[i].timeNs = (double)getTipFitTime(tip_hit,tip_waveform_pretrigger);
            sortedEvt.csiHit[i].PID = (float)tip_hit->GetPID();
            sortedEvt.csiHit[i].detNum = (uint8_t)tip_hit->GetTipChannel();
            fwrite(&sortedEvt.csiHit[i],sizeof(csi_hit),1,out);
          }

          //write footer value
          fwrite(&footerVal,sizeof(uint8_t),1,out);
          
          numSeparatedEvents++;
          
        }
      }

    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  //write the number of separated events to the beginning of the file
  numSepEvts = (uint64_t)numSeparatedEvents;
  fseek(out,0,SEEK_SET);
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  cout << "Number of separated events: " << numSepEvts << " ("<< 100.0f*numSepEvts/(float)analentries << "\% of total)" << endl;

  cout << "Wrote separated data to: " << outfile << endl;

  fclose(out);

  //verify the number of events actually written
  VerifyTimingSep(AnalysisTree,outfile,calfile);
  analysisfile->Close();

}
int main(int argc, char **argv){

  SeparatorTiming *mysort = new SeparatorTiming();

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting SeparatorTiming code\n");

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
    outfile = "Separated.smol";
  }else if (argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Separated.smol";
  }else if (argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
  }else{
    printf("Too many arguments\nArguments: SeparatorTiming analysis_tree calibration_file output_file\n");
    return 0;
  }

  printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
