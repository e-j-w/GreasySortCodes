//use the Makefile!

#define TTGamma_SMOL_cxx
#include "common.h"
#include "TTGamma_SMOL.h"

using namespace std;

int eGate[2], eGamma[2];

void TTGamma_SMOL::SortData(const char *sfile, const char *outfile)
{
  //setup ROOT hists
  TList *ttgList = new TList;
  TH1F *ttg = new TH1F("TIGRESS_Total_Rate", "TIGRESS Total Rate", 2048, -1024, 1024);
  ttg->GetXaxis()->SetTitle("t_{diff} (ns)");
  ttg->GetYaxis()->SetTitle("Counts / ns");
  ttgList->Add(ttg);

  //open SMOL file
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

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numTigHits; tigHitInd++){
      float addE1 = sortedEvt.tigHit[tigHitInd].energy;
      if(addE1 > MIN_TIG_EAB){
        if((addE1 >= eGate[0])&&(addE1 <= eGate[1])){
          for(int tigHitInd2 = 0; tigHitInd2 < sortedEvt.header.numTigHits; tigHitInd2++){
            if(tigHitInd2 != tigHitInd){
              float addE2 = sortedEvt.tigHit[tigHitInd2].energy;
              if(addE2 > MIN_TIG_EAB){
                if((addE2 >= eGamma[0])&&(addE2 <= eGamma[1])){
                  //cascade identified
                  Double_t tDiff = tigHitTime(&sortedEvt,tigHitInd) - tigHitTime(&sortedEvt,tigHitInd2);
                  /*if(tDiff < 0){
                    tDiff = -1.*tDiff;
                  }*/
                  ttg->Fill(tDiff);
                }
              }
            }
          }
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << endl << "Event sorting complete." << endl;
  fclose(inp);

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  ttgList->Write();
  myfile->Write();
  myfile->Close();

}
int main(int argc, char **argv)
{

  TTGamma_SMOL *mysort = new TTGamma_SMOL();

  const char *sfile;
  const char *outfile;
  printf("Starting TTGamma_SMOL\n");
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
    cout << "Code sorts gamma-gamma timing." << endl;
    cout << "Arguments: TTGamma_SMOL smol_file eLow1 eHigh1 eLow2 eHigh2 out_file" << endl;
    return 0;
  }else if(argc == 7){
    sfile = argv[1];
    eGate[0] = atoi(argv[2]);
    eGate[1] = atoi(argv[3]);
    eGamma[0] = atoi(argv[4]);
    eGamma[1] = atoi(argv[5]);
    outfile = argv[6];
    printf("SMOL file: %s\n", sfile); 
  }else{
    printf("ERROR: Improper number of arguments!\nArguments: TTGamma_SMOL smol_file eLow1 eHigh1 eLow2 eHigh2 out_file\n");
    return 0;
  }

  if(eGate[0] > eGate[1]){
    //swap values
    int swapVal = eGate[1];
    eGate[1] = eGate[0];
    eGate[0] = swapVal;
  }
  if(eGamma[0] > eGamma[1]){
    //swap values
    int swapVal = eGamma[1];
    eGamma[1] = eGamma[0];
    eGamma[0] = swapVal;
  }

  cout << "Gate 1: [" << eGate[0] << " " << eGate[1] << "] keV" << endl;
  cout << "Gate 2: [" << eGamma[0] << " " << eGamma[1] << "] keV" << endl;
  cout << "Output file: " << outfile << endl;

  mysort->SortData(sfile,outfile);

  return 0;
}
