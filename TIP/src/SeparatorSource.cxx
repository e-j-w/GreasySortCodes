//Sort code to separate source data (TIGRESS only)

#define SeparatorSource_cxx
#include "common.h"
#include "evt_fmt.h"
#include "SeparatorSource.h"

using namespace std;

FILE *out;

uint64_t SeparatorSource::SortData(const char *afile, const char *calfile, const int noAddback){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("WARNING: opening file %s failed!\n", afile);
    return 0;
  }

  printf("File %s opened\n", afile);
  TTree *AnalysisTree = (TTree *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if(AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "ERROR: Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
    return 0;
  }

  TTip *tip = 0;
  if(AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "No TIP data present." << endl;
  }

  TTigressHit *add_hit;
  TTipHit *tip_hit;
  sorted_evt sortedEvt;
  uint8_t footerVal = 227U;

  TChannel::ReadCalFile(calfile);

  uint64_t numSeparatedEvents = 0;

  printf("\nSorting analysis events...\n");
  for(uint64_t jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    if(tip && (tip->GetMultiplicity()>0)){
      continue;
    }

    if(tigress){

      if((noAddback==0)&&((tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT)||(tigress->GetAddbackMultiplicity()>MAX_EVT_HIT))){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
        continue;
      }else if((noAddback==1)&&((tigress->GetMultiplicity()>MAXNUMTIGHIT)||(tigress->GetMultiplicity()>MAX_EVT_HIT))){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetMultiplicity() << ")!" << endl;
        continue;
      }else{

        memset(&sortedEvt,0,sizeof(sorted_evt));
        
        uint8_t numTigHits = 0;
        if(noAddback==0){
          for(int i = 0; i<tigress->GetAddbackMultiplicity();i++){
            add_hit = tigress->GetAddbackHit(i);
            if(!(add_hit->BGOFired()) && (add_hit->GetEnergy() > 0)){
              if(sortedEvt.header.evtTimeNs == 0){
                sortedEvt.header.evtTimeNs = (double)add_hit->GetTime();
              }
              sortedEvt.tigHit[numTigHits].energy = (float)add_hit->GetEnergy();
              sortedEvt.tigHit[numTigHits].timeOffsetNs = (float)(add_hit->GetTime() - sortedEvt.header.evtTimeNs);
              sortedEvt.tigHit[numTigHits].core = (uint8_t)add_hit->GetArrayNumber();
              sortedEvt.tigHit[numTigHits].seg = (uint8_t)add_hit->GetFirstSeg();
              numTigHits++;
            }
          }
        }else{
          for(int i = 0; i<tigress->GetMultiplicity();i++){
            add_hit = tigress->GetTigressHit(i);
            if(!(add_hit->BGOFired()) && (add_hit->GetEnergy() > 0)){
              if(sortedEvt.header.evtTimeNs == 0){
                sortedEvt.header.evtTimeNs = (double)add_hit->GetTime();
              }
              sortedEvt.tigHit[numTigHits].energy = (float)add_hit->GetEnergy();
              sortedEvt.tigHit[numTigHits].timeOffsetNs = (float)(add_hit->GetTime() - sortedEvt.header.evtTimeNs);
              sortedEvt.tigHit[numTigHits].core = (uint8_t)add_hit->GetArrayNumber();
              sortedEvt.tigHit[numTigHits].seg = (uint8_t)add_hit->GetFirstSeg();
              numTigHits++;
            }
          }
        }

        sortedEvt.header.numTigHits = numTigHits;
        sortedEvt.header.numCsIHits = (uint8_t)0;
        sortedEvt.header.numRFHits = (uint8_t)0;
        fwrite(&sortedEvt.header,sizeof(evt_header),1,out);

        for(int i = 0; i<numTigHits;i++){
          fwrite(&sortedEvt.tigHit[i],sizeof(tig_hit),1,out);
        }
        //write footer value
        fwrite(&footerVal,sizeof(uint8_t),1,out);
        
        numSeparatedEvents++;

      }

    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  cout << "Number of separated events: " << numSeparatedEvents << " ("<< 100.0f*numSeparatedEvents/(float)analentries << "\% of total)" << endl;


  analysisfile->Close();

  return (uint64_t)numSeparatedEvents;

}
int main(int argc, char **argv){

  SeparatorSource *mysort = new SeparatorSource();

  char const *afile;
  char const *outfile;
  char const *calfile;
  int noAddback;
  printf("Starting SeparatorSource code\n");

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
    cout << "Arguments: SeparatorSource analysis_tree calibration_file output_file noAddback" << endl;
    cout << "*analysis_tree* can be a single analysis tree (extension .root), or a list of analysis trees (extension .list, one filepath per line)." << endl;
    cout << "Default values will be used if arguments (other than analysis_tree) are omitted." << endl;
    return 0;
  }else if (argc == 2){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Separated.smol";
    noAddback = 0;
  }else if (argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Separated.smol";
    noAddback = 0;
  }else if (argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    noAddback = 0;
  }else if (argc == 5){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    noAddback = atoi(argv[4]);
  }else{
    printf("Incorrect arguments\nArguments: SeparatorSource analysis_tree calibration_file output_file noAddback\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  if(noAddback){
    cout << "Will not sort addback energies." << endl;
  }else{
    cout << "Will sort addback energies." << endl;
  }

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get analysis tree or list file name." << endl;
    return 0;
  }

  //setup the output file
  out = fopen(outfile, "wb");
  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  if(strcmp(dot + 1, "root") == 0){
    printf("Analysis tree file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    numSepEvts += mysort->SortData(afile, calfile, noAddback);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("Analysis tree list: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(afile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << afile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, calfile, noAddback);
        }
      }
    }
  }else{
    cout << "ERROR: improper file extension for analysis tree or list (should be .root or .list)." << endl;
    return 0;
  }

  //write the number of separated events to the beginning of the file
  fseek(out,0,SEEK_SET);
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);
  cout << "Wrote " << numSepEvts << " separated events to: " << outfile << endl;
  fclose(out);

  return 0;
}
