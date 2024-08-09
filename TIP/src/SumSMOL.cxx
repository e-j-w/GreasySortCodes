//use the Makefile!

#define SumS_cxx
#include "common.h"
#include "SumSMOL.h"

using namespace std;

#define MAX_NUM_FILES 2

uint64_t Sum::SortData(const char *sfile, FILE *out)
{

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  printf("\nMerging events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    fwrite(&sortedEvt.header,sizeof(evt_header),1,out);
    for(int i = 0; i<sortedEvt.header.numTigHits;i++){
      fwrite(&sortedEvt.tigHit[i],sizeof(tig_hit),1,out);
    }
    for(int i = 0; i<sortedEvt.header.numNoABHits;i++){
      fwrite(&sortedEvt.noABHit[i],sizeof(tig_hit),1,out);
    }
    for(int i = 0; i<sortedEvt.header.numCsIHits;i++){
      fwrite(&sortedEvt.csiHit[i],sizeof(csi_hit),1,out);
    }
    for(int i = 0; i<numRFHits(&sortedEvt);i++){
      cout << "RF not implemented!!!" << endl;
    }
    //write footer value
    footerVal = 227U;
    fwrite(&footerVal,sizeof(uint8_t),1,out);

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  fclose(inp);

  return sentries;

}

int main(int argc, char **argv)
{

  Sum *mysort = new Sum();

  const char *sfile[MAX_NUM_FILES];
  const char *outfile;
  printf("Starting SumSMOL\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if(argc == 1){
    cout << "Code adds together 2 data files in the SMOL format" << endl;
    cout << "Arguments: SumSMOL smol_file1 smol_file2 output_file" << endl;
    cout << "Default values will be used if arguments (other than input SMOL files) are omitted." << endl;
    return 0;
  }else if(argc == 3){
    sfile[0] = argv[1];
    sfile[1] = argv[2];
    outfile = "out.smol";
    printf("Input files: %s, %s\nOutput file: %s\n", sfile[0], sfile[1], outfile); 
  }else if(argc == 4){
    sfile[0] = argv[1];
    sfile[1] = argv[2];
    outfile = argv[3];
    printf("Input files: %s, %s\nOutput file: %s\n", sfile[0], sfile[1], outfile); 
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: SumSMOL smol_file1 smol_file2 output_file\n");
    return 0;
  }

  uint8_t numFiles = 2;
  
  //setup the output file
  FILE *out = fopen(outfile, "wb");
  uint64_t numSepEvts = 0U;
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);

  for(uint8_t i=0;i<numFiles;i++){
    numSepEvts += mysort->SortData(sfile[i], out);
  }

  //write the number of separated events to the beginning of the file
  fseek(out,0,SEEK_SET);
  fwrite(&numSepEvts,sizeof(uint64_t),1,out);
  cout << "Wrote " << numSepEvts << " separated events to: " << outfile << endl;
  fclose(out);

  return 0;
}
