//Generates TIGRESS gamma ray spectra for PID and time separated data
//timing windows are defined in common.h
//PID gates in common.cxx

#define ETIP_mca_SMOL_cxx
#include "common.h"
#include "ETIP_mca_SMOL.h"

using namespace std;

void ETIP_mca_SMOL::WriteData(const char* outName){

  cout << "Writing gated histogram to: " << outName << endl;

  FILE *out;
  if((out = fopen(outName, "w")) == NULL){ //open the file
    cout << "ERROR: Cannot open the output file: " << outName << endl;
    return;
  }else{
    fwrite(&mcaOut,sizeof(mcaOut),1,out);
    fclose(out);
  }

}

void ETIP_mca_SMOL::SortData(char const *sfile, double MeVPerBin){

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numTigHits; tipHitInd++){

      if(sortedEvt.csiHit[tipHitInd].energy > 0.0){
        int eTIP = (int)(sortedEvt.csiHit[tipHitInd].energy/MeVPerBin);
        if(eTIP>=0 && eTIP<S32K){
          mcaOut[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)][eTIP]++;
          mcaOut[NTIPRING][eTIP]++;
        }
      }
      
    }

    if (jentry % 1000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;
  
  fclose(inp);
  
}

int main(int argc, char **argv){

  ETIP_mca_SMOL *mysort = new ETIP_mca_SMOL();

  const char *sfile;
  const char *outfile;
  double MeVPerBin = 1.0;
  printf("Starting ETIP_mca_SMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  if((argc != 3)&&(argc != 4)){
    cout << "Generates CsI ball mca spectra." << endl;
    cout << "Arguments: ETIP_mca_SMOL smol_file output_file MeV_per_bin" << endl;
    cout << "  *MeV_per_bin* defaults to 1 if not specified." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    if(argc > 3){
      MeVPerBin = atof(argv[3]);
    }
  }

  if(MeVPerBin <= 0.0){
    cout << "ERROR: Invalid keV/bin factor (" << MeVPerBin << ")!" << endl;
    return 0;
  }

  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  //single analysis tree
  cout << "SMOL tree: " << sfile << endl;
  cout << "Output file: " << outfile << endl;
  cout << "Written spectra will have " << MeVPerBin << " MeV per bin." << endl;

  mysort->SortData(sfile, MeVPerBin);
  
  mysort->WriteData(outfile);

  return 0;
}
