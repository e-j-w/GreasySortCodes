//use the Makefile!

#define SortDiagnosticsS_cxx
#include "common.cxx"
#include "SortDiagnosticsSMOL.h"

using namespace std;

sorted_evt sortedEvt;

void SortDiagnosticsS::SortData(const char *sfile, const char *outfile)
{
  Initialise();

  FILE *inp = fopen(sfile, "rb");
  if(inp==NULL){
    printf("ERROR: couldn't open file: %s\n",sfile);
    exit(-1);
  }
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  
  uint8_t footerVal;
  unsigned long int numHPGeABHits = 0;

  /*cout << "TIGRESS positions: " << endl;
  for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
      TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
      cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
  }*/

  printf("\nSorting events...\n");
  for(Long64_t jentry = 0; jentry < sentries; jentry++){

    //read event
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }

    hpgeMult->Fill(sortedEvt.header.numNoABHits);
    
    for (int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy);
        hpgeE_ANum->Fill(sortedEvt.noABHit[noABHitInd].core, sortedEvt.noABHit[noABHitInd].energy);
      }
    }

    //evaluate non-addback timing conditions
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){

        //TIG-TIG timing, position, and energy
        for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
          if(sortedEvt.noABHit[noABHitInd2].energy > MIN_HPGE_EAB){
            Double_t tDiff = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,noABHitInd2);
            hpgeT_hpgeT->Fill(tDiff);
            hpgeE_hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
            hpgeE_hpgeE->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
            hpgePos_hpgePos->Fill(sortedEvt.noABHit[noABHitInd].core,sortedEvt.noABHit[noABHitInd2].core);
            hpge_hpge_dist->Fill(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core,0,sortedEvt.noABHit[noABHitInd2].core,0,1)); //FORWARD POSITION (11 cm)
            hpge_hpge_angle->Fill(getGeVector(sortedEvt.noABHit[noABHitInd].core,0,1).Angle(getGeVector(sortedEvt.noABHit[noABHitInd2].core,0,1))*180.0/PI); //FORWARD POSITION (11 cm)
            if((tDiff >= hpgehpgeTGate[0])&&(tDiff <= hpgehpgeTGate[1])){
              hpgeT_hpgeT_tsep->Fill(tDiff);
              hpgeE_hpgeE_tsep->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
              hpgeE_hpgeE_tsep->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
              if(sortedEvt.header.numNoABHits == 2){
                hpgeT_hpgeT_tsepmult2->Fill(tDiff);
                hpgeE_hpgeE_tsepmult2->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
                hpgeE_hpgeE_tsepmult2->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
              }
            }
            if((tDiff >= hpgehpgeABGate[0])&&(tDiff <= hpgehpgeABGate[1])){
              if(getGeVector(sortedEvt.noABHit[noABHitInd].core,0,1).Angle(getGeVector(sortedEvt.noABHit[noABHitInd2].core,0,1))*180.0/PI > 175.0){
                hpgeE_hpgeE_180deg->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
                hpgeE_hpgeE_180deg->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
                hpgeE_hpgeE_180deg_proj->Fill(sortedEvt.noABHit[noABHitInd].energy);
                hpgeE_hpgeE_180deg_proj->Fill(sortedEvt.noABHit[noABHitInd2].energy);
                hpgeE_hpgeE_180deg_sum->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy);
              }
            }
            if((tDiff >= hpgehpgeTRandGate[0])&&(tDiff <= hpgehpgeTRandGate[1])){
              hpgeT_hpgeT_tseprand->Fill(tDiff);
              hpgeE_hpgeE_tseprand->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
              hpgeE_hpgeE_tseprand->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
            }
          }
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Number of HPGe addback hits: " << numHPGeABHits << endl;
  cout << endl << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *hpgedir = myfile->mkdir("HPGe");
  hpgedir->cd();
  hpgeList->Write();
  myfile->cd();

  TDirectory *timingdir = myfile->mkdir("Timing");
  timingdir->cd();
  timingList->Write();
  myfile->cd();

  TDirectory *hpgehpgedir = myfile->mkdir("HPGe_HPGe");
  hpgehpgedir->cd();
  hpgehpgeList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
  fclose(inp);
}
int main(int argc, char **argv)
{

  SortDiagnosticsS *mysort = new SortDiagnosticsS();

  const char *sfile;
  const char *outfile;
  printf("Starting SortDiagnosticsSMOL\n");

  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online HPGe data" << endl;
    cout << "Arguments: SortDiagnosticsSMOL smol_file output_file" << endl;
    cout << "Default values will be used if arguments (other than SMOL file) are omitted." << endl;
    return 0;
  }else if(argc == 2){
    sfile = argv[1];
    outfile = "Histograms.root";
    printf("SMOL file: %s\nOutput file: %s\n", sfile, outfile); 
  }else if(argc == 3){
    sfile = argv[1];
    outfile = argv[2];
    printf("SMOL file: %s\nOutput file: %s\n", sfile, outfile);
  }else{
    printf("ERROR: too many arguments!\nArguments: SortDiagnostics SMOL smol_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(sfile, outfile);

  return 0;
}
