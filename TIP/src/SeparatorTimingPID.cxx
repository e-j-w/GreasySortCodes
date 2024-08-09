//Sort code to separate data based on TIP and TIGRESS timing
//as well as PID
//Timing windows are defined in common.h

#define SeparatorTimingPID_cxx
#include "common.h"
#include "evt_fmt.h"
#include "SeparatorTimingPID.h"

using namespace std;

FILE *out;

uint64_t SeparatorTiming::SortData(const char *afile, const int nP, const int nA, const char *calfile, PIDGates *gates){

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
    cout << "ERROR: Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
    return 0;
  }

  TTigressHit *add_hit, *noAB_hit;
  TTipHit *tip_hit;
  sorted_evt sortedEvt;
  uint8_t footerVal = 227U;

  TChannel::ReadCalFile(calfile);

  uint64_t numSeparatedEvents = 0;
  uint8_t nPTimingCondition, nATimingCondition;
  if(nP < 0){
    nPTimingCondition = 0;
  }else{
    nPTimingCondition = (uint8_t)nP;
  }
  if(nA < 0){
    nATimingCondition = 0;
  }else{
    nATimingCondition = (uint8_t)nA;
  }

  printf("\nSorting analysis events...\n");
  for(uint64_t jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    if(tigress && tip){

      tigress->ResetAddback();

      if((tip->GetMultiplicity()>MAXNUMTIPHIT)||(tip->GetMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
        continue;
      }else if((tigress->GetMultiplicity()>MAXNUMTIGHIT)||(tigress->GetMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetMultiplicity() << ")!" << endl;
        continue;
      }else{

        uint64_t passedtimeGateNoAB = passesTimeGateNoAB(tigress,tip,1,nPTimingCondition+nATimingCondition); //also rejects pileup
        
        if(passedtimeGateNoAB&(1ULL<<TIPTIGFLAG)){

          uint8_t evtNumAlphas = 0;
          uint8_t evtNumProtons = 0;
          for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
            if(passedtimeGateNoAB&(1ULL<<tipHitInd)){
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

          if(((evtNumProtons == nP)||(nP < 0))&&((evtNumAlphas == nA)||(nA < 0))){

            memset(&sortedEvt,0,sizeof(sorted_evt));
            
            uint8_t numTigHits = 0;
            uint8_t numNoABHits = 0;
            uint8_t suppressorFired = 0;
            for(int i = 0; i<tigress->GetMultiplicity();i++){
              if(i<MAX_EVT_HIT){
                noAB_hit = tigress->GetTigressHit(i);
                if(passedtimeGateNoAB&(1ULL<<(i+MAXNUMTIPHIT))){
                  if(!(noAB_hit->BGOFired()) && (noAB_hit->GetEnergy() > 0)){
                    if(sortedEvt.header.evtTimeNs == 0){
                      sortedEvt.header.evtTimeNs = (double)noAB_hit->GetTime();
                    }
                    sortedEvt.noABHit[numNoABHits].energy = (float)noAB_hit->GetEnergy();
                    sortedEvt.noABHit[numNoABHits].timeOffsetNs = (float)(noAB_hit->GetTime() - sortedEvt.header.evtTimeNs);
                    sortedEvt.noABHit[numNoABHits].core = (uint8_t)noAB_hit->GetArrayNumber();
                    sortedEvt.noABHit[numNoABHits].seg = (uint8_t)noAB_hit->GetFirstSeg();
                    numNoABHits++;
                  }
                }
                if(noAB_hit->BGOFired()){
                  suppressorFired = 1;
                }
              }
            }
            uint64_t passedtimeGate = passesTimeGateAB(tigress,tip,1,nPTimingCondition+nATimingCondition); //also rejects pileup
            for(int i = 0; i<tigress->GetAddbackMultiplicity();i++){
              add_hit = tigress->GetAddbackHit(i);
              if(passedtimeGate&(1ULL<<(i+MAXNUMTIPHIT))){
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
                if(add_hit->BGOFired()){
                  suppressorFired = 1;
                }
              }
            }
            uint8_t numCsIHits = 0;
            for(int i = 0; i<tip->GetMultiplicity();i++){
              tip_hit = tip->GetTipHit(i);
              if(passedtimeGate&(1ULL<<i)){
                if(sortedEvt.header.evtTimeNs == 0){
                  sortedEvt.header.evtTimeNs = (double)add_hit->GetTime();
                }
                sortedEvt.csiHit[numCsIHits].energy = (float)tip_hit->GetEnergy();
                sortedEvt.csiHit[numCsIHits].timeOffsetNs = (float)(getTipFitTime(tip_hit,tip_waveform_pretrigger) - sortedEvt.header.evtTimeNs);
                sortedEvt.csiHit[numCsIHits].PID = (float)tip_hit->GetPID();
                sortedEvt.csiHit[numCsIHits].detNum = (uint8_t)tip_hit->GetTipChannel();
                sortedEvt.csiHit[numCsIHits].metadata |= ((uint8_t)(tip_hit->GetFitType()) & 7U); //store fit type
                numCsIHits++;
              }
            }

            sortedEvt.header.numTigHits = numTigHits;
            sortedEvt.header.numNoABHits = numNoABHits;
            sortedEvt.header.numCsIHits = numCsIHits;
            sortedEvt.header.metadata = (uint8_t)0;
            if(suppressorFired){
              sortedEvt.header.metadata |= (uint8_t)(1U << 1);
            }
            fwrite(&sortedEvt.header,sizeof(evt_header),1,out);

            for(int i = 0; i<numTigHits;i++){
              fwrite(&sortedEvt.tigHit[i],sizeof(tig_hit),1,out);
            }

            for(int i = 0; i<numNoABHits;i++){
              fwrite(&sortedEvt.noABHit[i],sizeof(tig_hit),1,out);
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
  int nP = 0;
  int nA = 0;

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Arguments: SeparatorTimingPID analysis_tree numP numA calibration_file output_file" << endl;
    cout << "*analysis_tree* can be a single analysis tree (extension .root), or a list of analysis trees (extension .list, one filepath per line)." << endl;
    cout << "nP = num protons required, nA = num alphas required.  Values of -1 indicate no condition on that particle type" << endl;
    cout << "if no_addback = 1, will sort non-addback energies for TIGRESS" << endl;
    cout << "Default values will be used if arguments (other than analysis tree, nP, nA) are omitted." << endl;
    return 0;
  }else if (argc == 3){
    afile = argv[1];
    nP = atoi(argv[2]);
    nA = atoi(argv[3]);
    calfile = "CalibrationFile.cal";
    outfile = "Separated.smol";
  }else if (argc == 4){
    afile = argv[1];
    nP = atoi(argv[2]);
    nA = atoi(argv[3]);
    calfile = argv[4];
    outfile = "Separated.smol";
  }else if (argc == 5){
    afile = argv[1];
    nP = atoi(argv[2]);
    nA = atoi(argv[3]);
    calfile = argv[4];
    outfile = argv[5];
  }else if (argc == 6){
    afile = argv[1];
    nP = atoi(argv[2]);
    nA = atoi(argv[3]);
    calfile = argv[4];
    outfile = argv[5];
  }else{
    printf("Incorrect arguments\nArguments: SeparatorTimingPID analysis_tree numP numA calibration_file output_file\n");
    return 0;
  }

  printf("Starting SeparatorTimingPID code\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

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
    printf("Analysis tree file: %s\nNumber of protons: %i\nNumber of alphas: %i\nCalibration file: %s\nOutput file: %s\n", afile, nP, nA, calfile, outfile);
    numSepEvts += mysort->SortData(afile, nP, nA, calfile, gates);
  }else if(strcmp(dot + 1, "list") == 0){
    printf("Analysis tree list: %s\nNumber of protons: %i\nNumber of alphas: %i\nCalibration file: %s\nOutput file: %s\n", afile, nP, nA, calfile, outfile);
    
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
