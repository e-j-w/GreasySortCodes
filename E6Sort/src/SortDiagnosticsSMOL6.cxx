//use the Makefile!

#define SortDiagnosticsS6_cxx
#include "common.cxx"
#include "SortDiagnosticsSMOL6.h"

using namespace std;

sorted_evt sortedEvt;

void SortDiagnosticsS6::SortData(const char *sfile, const char *outfile)
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
    if(readSMOL6Event(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
    }
    
    for (int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
        tigE->Fill(sortedEvt.noABHit[tigHitInd].energy);
        tigE_ANum->Fill(sortedEvt.noABHit[tigHitInd].core, sortedEvt.noABHit[tigHitInd].energy);
        if(sortedEvt.header.numNoABHits == 2){
          tigE_mult2->Fill(sortedEvt.noABHit[tigHitInd].energy);
        }
      }
    }

    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){

      //cout << "energy: " << sortedEvt.tigHit[tigHitIndAB].energy << ", array num: " << sortedEvt.tigHit[tigHitIndAB].core << ", address: " << sortedEvt.tigHit[tigHitIndAB]->GetAddress() << endl;

      numTigressABHits++;

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
        addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy);
        tigRate->Fill(tigHitTime(&sortedEvt,tigHitIndAB)/pow(10,9));
        tigNum_time->Fill(tigHitTime(&sortedEvt,tigHitIndAB)/pow(10,9),sortedEvt.tigHit[tigHitIndAB].core);
        //cout << "hit " << tigHitIndAB << ", Anum: " << sortedEvt.tigHit[tigHitIndAB].core << ", energy: " << sortedEvt.tigHit[tigHitIndAB].energy << endl;
        addE_ANum->Fill(sortedEvt.tigHit[tigHitIndAB].core, sortedEvt.tigHit[tigHitIndAB].energy);
        TVector3 tigVec = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg);
        //tigVec.SetZ(tigVec.Z() - 5.0); //target position offset
        double theta = tigVec.Theta()*180./PI;
        double phi = tigVec.Phi()*180./PI;
        addE_theta->Fill(theta, sortedEvt.tigHit[tigHitIndAB].energy);
        addE_phi->Fill(phi, sortedEvt.tigHit[tigHitIndAB].energy);
        theta_phi->Fill(theta, phi);
        //fill ring spectra
        addE_ring[getTIGRESSRing(theta)]->Fill(sortedEvt.tigHit[tigHitIndAB].energy);
      }
      
    }

    //evaluate non-addback timing conditions
    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){

      if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){

        //TIG-TIG timing, position, and energy
        for(int tigHitInd2 = tigHitInd+1; tigHitInd2 < sortedEvt.header.numNoABHits; tigHitInd2++){
          if(sortedEvt.noABHit[tigHitInd2].energy > MIN_TIG_EAB){
            Double_t tDiff = noABHitTime(&sortedEvt,tigHitInd) - noABHitTime(&sortedEvt,tigHitInd2);
            tigT_tigT->Fill(tDiff);
            tigE_tigE->Fill(sortedEvt.noABHit[tigHitInd].energy,sortedEvt.noABHit[tigHitInd2].energy);
            tigE_tigE->Fill(sortedEvt.noABHit[tigHitInd2].energy,sortedEvt.noABHit[tigHitInd].energy); //symmetrized
            if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
              tigT_tigT_tsep->Fill(tDiff);
              tigE_tigE_tsep->Fill(sortedEvt.noABHit[tigHitInd].energy,sortedEvt.noABHit[tigHitInd2].energy);
              tigE_tigE_tsep->Fill(sortedEvt.noABHit[tigHitInd2].energy,sortedEvt.noABHit[tigHitInd].energy); //symmetrized
              if(sortedEvt.header.numNoABHits == 2){
                tigT_tigT_tsepmult2->Fill(tDiff);
                tigE_tigE_tsepmult2->Fill(sortedEvt.noABHit[tigHitInd].energy,sortedEvt.noABHit[tigHitInd2].energy);
                tigE_tigE_tsepmult2->Fill(sortedEvt.noABHit[tigHitInd2].energy,sortedEvt.noABHit[tigHitInd].energy); //symmetrized
              }
            }else if((tDiff >= tigtigTBGGate[0])&&(tDiff <= tigtigTBGGate[1])){
              tigE_tigE_trand->Fill(sortedEvt.noABHit[tigHitInd].energy,sortedEvt.noABHit[tigHitInd2].energy);
              tigE_tigE_trand->Fill(sortedEvt.noABHit[tigHitInd2].energy,sortedEvt.noABHit[tigHitInd].energy); //symmetrized
            }
          }
        }
      }
    }

    //evaluate addback timing conditions
    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){

        //TIG-TIG timing, position, and energy
        for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
          if(sortedEvt.tigHit[tigHitIndAB2].energy > MIN_TIG_EAB){
            Double_t tDiff = tigHitTime(&sortedEvt,tigHitIndAB) - tigHitTime(&sortedEvt,tigHitIndAB2);
            addT_addT->Fill(tDiff);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
            tigPos_tigPos->Fill(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB2].core);
            if((tDiff >= tigtigTGate[0])&&(tDiff <= tigtigTGate[1])){
              addT_addT_tsep->Fill(tDiff);
              addE_addE_tsep->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
              addE_addE_tsep->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
              if(sortedEvt.header.numTigHits == 2){
                addT_addT_tsepmult2->Fill(tDiff);
                addE_addE_tsepmult2->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
                addE_addE_tsepmult2->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
              }
            }else if((tDiff >= tigtigTBGGate[0])&&(tDiff <= tigtigTBGGate[1])){
              addE_addE_trand->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
              addE_addE_trand->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
            }
          }
        }
      }
    }

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigressABHits << endl;
  cout << endl << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  myfile->cd();

  TDirectory *timingdir = myfile->mkdir("Timing");
  timingdir->cd();
  timingList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS_TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
  fclose(inp);
}
int main(int argc, char **argv)
{

  SortDiagnosticsS6 *mysort = new SortDiagnosticsS6();

  const char *sfile;
  const char *outfile;
  printf("Starting SortDiagnosticsSMOL6\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online source data." << endl;
    cout << "Arguments: SortDiagnosticsSMOL6 smol6_file output_file" << endl;
    cout << "Default values will be used if arguments (other than SMOL6 file) are omitted." << endl;
    return 0;
  }else if(argc == 2){
    sfile = argv[1];
    outfile = "Histograms.root";
    printf("SMOL6 file: %s\nOutput file: %s\n", sfile, outfile); 
  }else if(argc == 3){
    sfile = argv[1];
    outfile = argv[2];
    printf("SMOL6 file: %s\nOutput file: %s\n", sfile, outfile);
  }else{
    printf("ERROR: too many arguments!\nExpected arguments: SortDiagnosticsSMOL6 smol6_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(sfile, outfile);

  return 0;
}
