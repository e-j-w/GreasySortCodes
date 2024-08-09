//use the Makefile!

#define TIPRingPairs_cxx
#include "common.h"
#include "TIPRingPairs.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];

void TIPRingPairs::SortData(const char *sfile)
{

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  unsigned long int numTigressABHits = 0;
  unsigned long int numTipHits = 0;
  unsigned long int numTipTigHits = 0;

  uint64_t tipRingCts[NTIPRING][NTIPRING];
  memset(tipRingCts,0,sizeof(tipRingCts));

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    //TIP-TIP timing, position, and energy
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
      for(int tipHitInd2 = tipHitInd+1; tipHitInd2 < sortedEvt.header.numCsIHits; tipHitInd2++){
        tipRingCts[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)][getTIPRing(sortedEvt.csiHit[tipHitInd2].detNum)]++;
        tipRingCts[getTIPRing(sortedEvt.csiHit[tipHitInd2].detNum)][getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)]++;
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigressABHits << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP + TIGRESS hits: " << numTipTigHits << endl << endl;
  cout << endl << "Event sorting complete" << endl;

  cout << "Ring 1   Ring 2   Counts" << endl;
  for(int i=0; i<NTIPRING; i++){
    for(int j=0; j<NTIPRING; j++){
      cout << i << " " << j << " " << tipRingCts[i][j] << endl;
    }
  }

  fclose(inp);
}
int main(int argc, char **argv)
{

  TIPRingPairs *mysort = new TIPRingPairs();

  const char *sfile;
  printf("Starting TIPRingPairs\n");
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0){
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Computes the number of hits in each pair of rings of the of TIP CsI ball." << endl;
    cout << "Arguments: TIPRingPairs smol_file" << endl;
    return 0;
  }else if(argc == 2){
    sfile = argv[1];
    printf("SMOL file: %s\n\n", sfile); 
  }else{
    printf("ERROR: wrong number of arguments!\nArguments: TIPRingPairs SMOL smol_file\n");
    return 0;
  }

  mysort->SortData(sfile);

  return 0;
}
