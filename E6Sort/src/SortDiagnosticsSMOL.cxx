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
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  
  uint8_t footerVal;
  unsigned long int numTigressABHits = 0;

  Int_t evtNumProtons, evtNumAlphas;
  Int_t evtNumProtonsDetSumGate, evtNumAlphasDetSumGate;
  Int_t evtNumHorDetSumGate;

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

    //correct GRIFFIN cores
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].core >= 64){
        sortedEvt.noABHit[noABHitInd].core = (uint8_t)(sortedEvt.noABHit[noABHitInd].core - 64);
      }
    }
    for(int ABHitInd = 0; ABHitInd < sortedEvt.header.numABHits; ABHitInd++){
      if(sortedEvt.ABHit[ABHitInd].core >= 64){
        sortedEvt.ABHit[ABHitInd].core = (uint8_t)(sortedEvt.ABHit[ABHitInd].core - 64);
      }
    }
    
    for (int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy);
        hpgeE_ANum->Fill(sortedEvt.noABHit[noABHitInd].core, sortedEvt.noABHit[noABHitInd].energy);
      }
    }

    for(int ABHitInd = 0; ABHitInd < sortedEvt.header.numABHits; ABHitInd++){

      //cout << "energy: " << sortedEvt.ABHit[ABHitInd].energy << ", array num: " << sortedEvt.ABHit[ABHitInd].core << ", address: " << sortedEvt.ABHit[ABHitInd]->GetAddress() << endl;

      numTigressABHits++;

      if(sortedEvt.ABHit[ABHitInd].energy > MIN_HPGE_EAB){
        addE->Fill(sortedEvt.ABHit[ABHitInd].energy);
        hpgeRate->Fill(ABHitTime(&sortedEvt,ABHitInd)/pow(10,9));
        hpgeNum_time->Fill(ABHitTime(&sortedEvt,ABHitInd)/pow(10,9),sortedEvt.ABHit[ABHitInd].core);
        //cout << "hit " << ABHitInd << ", Anum: " << sortedEvt.ABHit[ABHitInd].core << ", energy: " << sortedEvt.ABHit[ABHitInd].energy << endl;
        addE_ANum->Fill(sortedEvt.ABHit[ABHitInd].core, sortedEvt.ABHit[ABHitInd].energy);
        TVector3 tigVec = getTigVector(sortedEvt.ABHit[ABHitInd].core,sortedEvt.ABHit[ABHitInd].seg);
        //tigVec.SetZ(tigVec.Z() - 5.0); //target position offset
        double theta = tigVec.Theta()*180./PI;
        double phi = tigVec.Phi()*180./PI;
        addE_theta->Fill(theta, sortedEvt.ABHit[ABHitInd].energy);
        addE_phi->Fill(phi, sortedEvt.ABHit[ABHitInd].energy);
        theta_phi->Fill(theta, phi);
        //fill ring spectra
        addE_ring[getTIGRESSRing(theta)]->Fill(sortedEvt.ABHit[ABHitInd].energy);
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
          }
        }
      }
    }

    //evaluate addback timing conditions
    for(int ABHitInd = 0; ABHitInd < sortedEvt.header.numABHits; ABHitInd++){

      if(sortedEvt.ABHit[ABHitInd].energy > MIN_HPGE_EAB){

        //TIG-TIG timing, position, and energy
        for(int ABHitInd2 = ABHitInd+1; ABHitInd2 < sortedEvt.header.numABHits; ABHitInd2++){
          if(sortedEvt.ABHit[ABHitInd2].energy > MIN_HPGE_EAB){
            Double_t tDiff = ABHitTime(&sortedEvt,ABHitInd) - ABHitTime(&sortedEvt,ABHitInd2);
            addT_addT->Fill(tDiff);
            addE_addE->Fill(sortedEvt.ABHit[ABHitInd].energy,sortedEvt.ABHit[ABHitInd2].energy);
            addE_addE->Fill(sortedEvt.ABHit[ABHitInd2].energy,sortedEvt.ABHit[ABHitInd].energy); //symmetrized
            hpgePos_hpgePos->Fill(sortedEvt.ABHit[ABHitInd].core,sortedEvt.ABHit[ABHitInd2].core);
            if((tDiff >= hpgehpgeTGate[0])&&(tDiff <= hpgehpgeTGate[1])){
              addT_addT_tsep->Fill(tDiff);
              addE_addE_tsep->Fill(sortedEvt.ABHit[ABHitInd].energy,sortedEvt.ABHit[ABHitInd2].energy);
              addE_addE_tsep->Fill(sortedEvt.ABHit[ABHitInd2].energy,sortedEvt.ABHit[ABHitInd].energy); //symmetrized
              if(sortedEvt.header.numABHits == 2){
                addT_addT_tsepmult2->Fill(tDiff);
                addE_addE_tsepmult2->Fill(sortedEvt.ABHit[ABHitInd].energy,sortedEvt.ABHit[ABHitInd2].energy);
                addE_addE_tsepmult2->Fill(sortedEvt.ABHit[ABHitInd2].energy,sortedEvt.ABHit[ABHitInd].energy); //symmetrized
              }
            }
          }
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Number of HPGe addback hits: " << numTigressABHits << endl;
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

  // Input-chain-file, output-histogram-file
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