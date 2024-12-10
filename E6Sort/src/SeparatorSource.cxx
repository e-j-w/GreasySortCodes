//Sort code to separate source data (TIGRESS only)

#define SeparatorSource_cxx
#include "common.h"
#include "evt_fmt.h"
#include "SeparatorSource.h"

using namespace std;

FILE *out;

uint64_t SeparatorSource::SortData(const char *afile, const char *calfile){

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("WARNING: opening file %s failed!\n", afile);
    return 0;
  }

  printf("File %s opened\n", afile);
  TTree *AnalysisTree = (TTree *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  TGriffin *griffin = 0;
  TGriffinBgo *griffin_bgo = 0;
  if(AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
    cout << "Sorting TIGRESS data." << endl;
  }else if(AnalysisTree->FindBranch("TGriffin")){
    AnalysisTree->SetBranchAddress("TGriffin", &griffin);
    cout << "Sorting GRIFFIN data." << endl;
  }else{
    cout << "ERROR: Neither 'TTigress' or 'TGriffin' branches were found!" << endl;
    return 0;
  }

  if(AnalysisTree->FindBranch("TGriffinBgo")){
    AnalysisTree->SetBranchAddress("TGriffinBgo", &griffin_bgo);
    cout << "GRIFFIN BGO data is present." << endl;
  }

  TTigressHit *add_hit, *noAB_hit;
  TGriffinHit *add_hit_grif, *noAB_hit_grif;
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

    if(tigress){

      tigress->ResetAddback();

      if((tigress->GetMultiplicity()>MAXNUMHPGEHIT)||(tigress->GetMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many TIGRESS hits (" << tigress->GetMultiplicity() << ")!" << endl;
        continue;
      }else{

        memset(&sortedEvt,0,sizeof(sorted_evt));
        
        uint8_t numABHits = 0;
        uint8_t numNoABHits = 0;
        uint8_t suppressorFired = 0;
        for(int i = 0; i<tigress->GetMultiplicity();i++){
          if(i<MAX_EVT_HIT){
            noAB_hit = tigress->GetTigressHit(i);
            if(!(noAB_hit->BGOFired()) && (noAB_hit->GetEnergy() > 0)){
              if(sortedEvt.header.evtTimeNs == 0){
                sortedEvt.header.evtTimeNs = (double)noAB_hit->GetTime();
              }
              sortedEvt.noABHit[numNoABHits].energy = (float)noAB_hit->GetEnergy();
              sortedEvt.noABHit[numNoABHits].timeOffsetNs = (float)(noAB_hit->GetTime() - sortedEvt.header.evtTimeNs);
              sortedEvt.noABHit[numNoABHits].core = (uint8_t)noAB_hit->GetArrayNumber();
              if(sortedEvt.ABHit[numNoABHits].core >= 64){
                cout << "WARNING: invalid TIGRESS core: " << sortedEvt.ABHit[numABHits].core << endl;
                continue;
              }
              sortedEvt.noABHit[numNoABHits].seg = (uint8_t)noAB_hit->GetFirstSeg();
              numNoABHits++;
            }
            if(noAB_hit->BGOFired()){
              suppressorFired = 1;
            }
          }
        }
        for(int i = 0; i<tigress->GetAddbackMultiplicity();i++){
          add_hit = tigress->GetAddbackHit(i);
          if(!(add_hit->BGOFired()) && (add_hit->GetEnergy() > 0)){
            if(sortedEvt.header.evtTimeNs == 0){
              sortedEvt.header.evtTimeNs = (double)add_hit->GetTime();
            }
            sortedEvt.ABHit[numABHits].energy = (float)add_hit->GetEnergy();
            sortedEvt.ABHit[numABHits].timeOffsetNs = (float)(add_hit->GetTime() - sortedEvt.header.evtTimeNs);
            sortedEvt.ABHit[numABHits].core = (uint8_t)add_hit->GetArrayNumber();
            if(sortedEvt.ABHit[numABHits].core >= 64){
              cout << "WARNING: invalid TIGRESS core: " << sortedEvt.ABHit[numABHits].core << endl;
              continue;
            }
            sortedEvt.ABHit[numABHits].seg = (uint8_t)add_hit->GetFirstSeg();
            numABHits++;
          }
          if(add_hit->BGOFired()){
            suppressorFired = 1;
          }
        }

        sortedEvt.header.numABHits = numABHits;
        sortedEvt.header.numNoABHits = numNoABHits;
        sortedEvt.header.metadata = (uint8_t)0;
        if(suppressorFired){
          sortedEvt.header.metadata |= (uint8_t)(1U << 0);
        }
        fwrite(&sortedEvt.header,sizeof(evt_header),1,out);

        for(int i = 0; i<numABHits;i++){
          fwrite(&sortedEvt.ABHit[i].timeOffsetNs,sizeof(float),1,out);
          fwrite(&sortedEvt.ABHit[i].energy,sizeof(float),1,out);
          fwrite(&sortedEvt.ABHit[i].core,sizeof(uint8_t),1,out);
          fwrite(&sortedEvt.ABHit[i].seg,sizeof(uint8_t),1,out);
        }
        for(int i = 0; i<numNoABHits;i++){
          fwrite(&sortedEvt.noABHit[i].timeOffsetNs,sizeof(float),1,out);
          fwrite(&sortedEvt.noABHit[i].energy,sizeof(float),1,out);
          fwrite(&sortedEvt.noABHit[i].core,sizeof(uint8_t),1,out);
          fwrite(&sortedEvt.noABHit[i].seg,sizeof(uint8_t),1,out);
        }
        //write footer value
        fwrite(&footerVal,sizeof(uint8_t),1,out);
        
        numSeparatedEvents++;

      }

    }else if(griffin && griffin_bgo){

      griffin->ResetAddback();

      if((griffin->GetMultiplicity()>MAXNUMHPGEHIT)||(griffin->GetMultiplicity()>MAX_EVT_HIT)){
        cout << "WARNING: event " << jentry << " has too many GRIFFIN hits (" << griffin->GetMultiplicity() << ")!" << endl;
        continue;
      }else{

        memset(&sortedEvt,0,sizeof(sorted_evt));
        
        uint8_t numABHits = 0;
        uint8_t numNoABHits = 0;
        uint8_t suppressorFired = 0;
        uint64_t suppressorHitmap = 0;

        for(int i = 0; i<griffin->GetSuppressedMultiplicity(griffin_bgo);i++){
          if(i<MAX_EVT_HIT){
            noAB_hit_grif = griffin->GetSuppressedHit(i);
            if(noAB_hit_grif->GetEnergy() > 0){
              if(sortedEvt.header.evtTimeNs == 0){
                sortedEvt.header.evtTimeNs = (double)noAB_hit_grif->GetTime();
              }
              sortedEvt.noABHit[numNoABHits].energy = (float)noAB_hit_grif->GetEnergy();
              sortedEvt.noABHit[numNoABHits].timeOffsetNs = (float)(noAB_hit_grif->GetTime() - sortedEvt.header.evtTimeNs);
              sortedEvt.noABHit[numNoABHits].core = (uint8_t)(noAB_hit_grif->GetArrayNumber() + 63); //offset to indicate that this hit is GRIFFIN rather than TIGRESS, taking into account that GRIFFIN array number is 1-indexed while TIGRESS is 0-indexed
              numNoABHits++;
            }
          }
        }
        for(int i = 0; i<griffin->GetSuppressedAddbackMultiplicity(griffin_bgo);i++){
          add_hit_grif = griffin->GetSuppressedAddbackHit(i);
          if(add_hit_grif->GetEnergy() > 0){
            if(sortedEvt.header.evtTimeNs == 0){
              sortedEvt.header.evtTimeNs = (double)add_hit_grif->GetTime();
            }
            sortedEvt.ABHit[numABHits].energy = (float)add_hit_grif->GetEnergy();
            sortedEvt.ABHit[numABHits].timeOffsetNs = (float)(add_hit_grif->GetTime() - sortedEvt.header.evtTimeNs);
            sortedEvt.ABHit[numABHits].core = (uint8_t)(add_hit_grif->GetArrayNumber() + 63); //offset to indicate that this hit is GRIFFIN rather than TIGRESS, taking into account that GRIFFIN array number is 1-indexed while TIGRESS is 0-indexed
            numABHits++;
          }
        }

        sortedEvt.header.numABHits = numABHits;
        sortedEvt.header.numNoABHits = numNoABHits;
        sortedEvt.header.metadata = (uint8_t)0;
        if(griffin_bgo->GetMultiplicity() > 0){
          //at least one suppressor fired
          sortedEvt.header.metadata |= (uint8_t)(1U << 0);
        }
        fwrite(&sortedEvt.header,sizeof(evt_header),1,out);
        //write hits, without segment data (no segments for GRIFFIN)
        for(int i = 0; i<numABHits;i++){
          fwrite(&sortedEvt.ABHit[i].timeOffsetNs,sizeof(float),1,out);
          fwrite(&sortedEvt.ABHit[i].energy,sizeof(float),1,out);
          fwrite(&sortedEvt.ABHit[i].core,sizeof(uint8_t),1,out);
        }
        for(int i = 0; i<numNoABHits;i++){
          fwrite(&sortedEvt.noABHit[i].timeOffsetNs,sizeof(float),1,out);
          fwrite(&sortedEvt.noABHit[i].energy,sizeof(float),1,out);
          fwrite(&sortedEvt.noABHit[i].core,sizeof(uint8_t),1,out);
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

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Arguments: SeparatorSource analysis_tree calibration_file output_file" << endl;
    cout << "*analysis_tree* can be a single analysis tree (extension .root), or a list of analysis trees (extension .list, one filepath per line)." << endl;
    cout << "Default values will be used if arguments (other than analysis_tree) are omitted." << endl;
    return 0;
  }else if (argc == 1){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Separated.smole6";
  }else if (argc == 2){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Separated.smole6";
  }else if (argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
  }else if (argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
  }else{
    printf("Incorrect arguments\nArguments: SeparatorSource analysis_tree calibration_file output_file\n");
    return 0;
  }

  printf("Starting SeparatorSource code\n");

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

  //setup the output file
  out = fopen(outfile, "wb");
  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  if(strcmp(dot + 1, "root") == 0){
    printf("Analysis tree file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    numSepEvts += mysort->SortData(afile, calfile);
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
          numSepEvts += mysort->SortData(str, calfile);
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
