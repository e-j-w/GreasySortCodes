//Sort code to separate data based on TIP and TIGRESS timing
//as well as PID
//Timing windows are defined in common.h

#define SeparatorTimingPID_cxx
#include "common.h"
#include "evt_fmt.h"
#include "SeparatorTimingPID.h"

using namespace std;

uint64_t *mapping;
FILE *out;

uint64_t SeparatorTiming::SortData(const char *afile, const uint8_t nP, const uint8_t nA, const char *calfile, PIDGates *gates){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("WARNING: opening file %s failed!\n", afile);
    return 0;
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
    return 0;
  }

  TTip *tip = 0;
  if(AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "ERROR: Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
    return 0;
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

          uint8_t evtNumAlphas = 0;
          uint8_t evtNumProtons = 0;
          for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
            if(passedtimeGate&(1ULL<<tipHitInd)){
              tip_hit = tip->GetTipHit(tipHitInd);

              switch(getParticleType(tip_hit,gates)){ //see common.cxx
                case 4:
                  evtNumAlphas++;
                  break;
                case 1:
                  evtNumProtons++;
                  break;
                case 0:
                default:
                  break;
              }
            }
          }

          if((evtNumProtons == nP)&&(evtNumAlphas == nA)){
            
            mapping[numSeparatedEvents] = jentry;
            memset(&sortedEvt,0,sizeof(sorted_evt));
            
            uint8_t numTigHits = 0;
            for(int i = 0; i<tigress->GetAddbackMultiplicity();i++){
              add_hit = tigress->GetAddbackHit(i);
              if(passedtimeGate&(1ULL<<(i+MAXNUMTIPHIT))){
                if(!(add_hit->BGOFired()) && (add_hit->GetEnergy() > 0)){
                  sortedEvt.tigHit[numTigHits].energy = (float)add_hit->GetEnergy();
                  sortedEvt.tigHit[numTigHits].timeNs = (double)add_hit->GetTime();
                  sortedEvt.tigHit[numTigHits].core = (uint8_t)add_hit->GetArrayNumber();
                  sortedEvt.tigHit[numTigHits].seg = (uint8_t)add_hit->GetFirstSeg();
                  numTigHits++;
                }
              }
            }

            uint8_t numCsIHits = 0;
            for(int i = 0; i<tip->GetMultiplicity();i++){
              tip_hit = tip->GetTipHit(i);
              if(passedtimeGate&(1ULL<<i)){
                sortedEvt.csiHit[numCsIHits].energy = (float)tip_hit->GetEnergy();
                sortedEvt.csiHit[numCsIHits].timeNs = (double)getTipFitTime(tip_hit,tip_waveform_pretrigger);
                sortedEvt.csiHit[numCsIHits].PID = (float)tip_hit->GetPID();
                sortedEvt.csiHit[numCsIHits].detNum = (uint8_t)tip_hit->GetTipChannel();
                numCsIHits++;
              }
            }

            sortedEvt.header.numTigABHits = numTigHits;
            sortedEvt.header.numCsIHits = numCsIHits;
            sortedEvt.header.numRFHits = (uint8_t)0;
            fwrite(&sortedEvt.header,sizeof(evt_header),1,out);

            for(int i = 0; i<numTigHits;i++){
              fwrite(&sortedEvt.tigHit[i],sizeof(tigab_hit),1,out);
            }

            for(int i = 0; i<numCsIHits;i++){
              fwrite(&sortedEvt.csiHit[i],sizeof(csi_hit),1,out);
            }

            //write footer value
            fwrite(&footerVal,sizeof(uint8_t),1,out);
            
            numSeparatedEvents++;
          }
          
        }
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

  SeparatorTiming *mysort = new SeparatorTiming();

  char const *afile;
  char const *outfile;
  char const *calfile;
  uint8_t nP = 0;
  uint8_t nA = 2;
  printf("Starting SeparatorTimingPID code\n");

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
    cout << "Arguments: SeparatorTimingPID analysis_tree numP numA calibration_file output_file" << endl;
    cout << "*analysis_tree* can be a single analysis tree (extension .root), or a list of analysis trees (extension .list, one filepath per line)." << endl;
    cout << "Default values will be used if arguments (other than analysis tree, nP, nA) are omitted." << endl;
    return 0;
  }else if (argc == 4){
    afile = argv[1];
    nP = (uint8_t)atoi(argv[2]);
    nA = (uint8_t)atoi(argv[3]);
    calfile = "CalibrationFile.cal";
    outfile = "Separated.smol";
  }else if (argc == 5){
    afile = argv[1];
    nP = (uint8_t)atoi(argv[2]);
    nA = (uint8_t)atoi(argv[3]);
    calfile = argv[4];
    outfile = "Separated.smol";
  }else if (argc == 6){
    afile = argv[1];
    nP = (uint8_t)atoi(argv[2]);
    nA = (uint8_t)atoi(argv[3]);
    calfile = argv[4];
    outfile = argv[5];
  }else{
    printf("Incorrect arguments\nArguments: SeparatorTimingPID analysis_tree numP numA calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  const char *dot = strrchr(afile, '.'); //get the file extension
  if(dot==NULL){
    cout << "ERROR: couldn't get analysis tree or list file name." << endl;
    return 0;
  }

  //Setup TIP PID gates
  PIDGates *gates = new PIDGates;

  //setup the output file
  out = fopen(outfile, "wb");
  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  if(strcmp(dot + 1, "root") == 0){
    printf("Analysis tree file: %s\nNumber of protons: %u\nNumber of alphas: %u\nCalibration file: %s\nOutput file: %s\n", afile, nP, nA, calfile, outfile);
    numSepEvts += mysort->SortData(afile, nP, nA, calfile, gates);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("Analysis tree list: %s\nNumber of protons: %u\nNumber of alphas: %u\nCalibration file: %s\nOutput file: %s\n", afile, nP, nA, calfile, outfile);
    
    FILE *listfile;
    char str[256];

    if((listfile=fopen(afile,"r"))==NULL){
      cout << "ERROR: Cannot open the list file: " << afile << endl;
      return 0;
    }else{
      while(!(feof(listfile))){//go until the end of file is reached
        if(fgets(str,256,listfile)!=NULL){ //get an entire line
          str[strcspn(str, "\r\n")] = 0;//strips newline characters from the string
          numSepEvts += mysort->SortData(str, nP, nA, calfile, gates);
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
