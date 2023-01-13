//use the Makefile!

#define SortDiagnosticsS_cxx
#include "common.h"
#include "SortDiagnosticsSMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];

void SortDiagnostics::SortData(char const *sfile, char const *outfile)
{
  Initialise();

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  sorted_evt sortedEvt;
  uint8_t footerVal;

  unsigned long int numTigressABHits = 0;
  unsigned long int numTipHits = 0;
  unsigned long int numTipTigHits = 0;

  double tipFitTimes[MAXNUMTIPHIT];

  double tipEEvt;
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
    memset(&sortedEvt,0,sizeof(sorted_evt));
    footerVal = 0;
    fread(&sortedEvt.header,sizeof(evt_header),1,inp);
    for(int i = 0; i<sortedEvt.header.numTigABHits;i++){
      fread(&sortedEvt.tigHit[i],sizeof(tigab_hit),1,inp);
    }
    for(int i = 0; i<sortedEvt.header.numCsIHits;i++){
      fread(&sortedEvt.csiHit[i],sizeof(csi_hit),1,inp);
    }
    fread(&footerVal,sizeof(uint8_t),1,inp);
    if(footerVal != 227U){
      printf("ERROR: invalid footer value in event %lu (%u)!\n", jentry, footerVal);
      exit(-1);
    }

    if(sortedEvt.header.numCsIHits>MAXNUMTIPHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIP hits (" << sortedEvt.header.numCsIHits << ")!" << endl;
      continue;
    }
    if(sortedEvt.header.numTigABHits>MAXNUMTIGHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIGRESS hits (" << sortedEvt.header.numTigABHits << ")!" << endl;
      continue;
    }
    
    evtNumProtons=0;
    evtNumAlphas=0;
    evtNumProtonsDetSumGate=0;
    evtNumAlphasDetSumGate=0;
    evtNumHorDetSumGate=0;
    
      
    if(sortedEvt.header.numCsIHits>=0){
      tip_mult->Fill(sortedEvt.header.numCsIHits);
    }

    tipEEvt = 0.;
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){

      if((sortedEvt.csiHit[tipHitInd].detNum > 0)&&(sortedEvt.csiHit[tipHitInd].detNum <= NTIP)){

        tipRate->Fill(sortedEvt.csiHit[tipHitInd].timeNs/pow(10,9));
        
        numTipRingHits[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)]++;
        numTipHits++;
        
        tip_E->Fill(sortedEvt.csiHit[tipHitInd].energy);
        tipEEvt += sortedEvt.csiHit[tipHitInd].energy;
        tip_pos->Fill(sortedEvt.csiHit[tipHitInd].detNum);
        tip_ring->Fill(getTIPRing(sortedEvt.csiHit[tipHitInd].detNum));
        tip_E_pos->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.csiHit[tipHitInd].detNum);

        //PID related stuff
        if(sortedEvt.csiHit[tipHitInd].PID >= 0){ //PID was found
          tip_E_PID_Sum->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.csiHit[tipHitInd].PID);
          //cout << "Energy: " << sortedEvt.csiHit[tipHitInd].energy << ", PID: " << sortedEvt.csiHit[tipHitInd].PID << endl;
          tip_E_PID[sortedEvt.csiHit[tipHitInd].detNum-1]->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.csiHit[tipHitInd].PID);
          tip_E_PID_Ring[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)]->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.csiHit[tipHitInd].PID);
        }
      }
      
    }
    tip_Etot->Fill(tipEEvt);

    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigABHits; tigHitIndAB++){

      //cout << "energy: " << sortedEvt.tigHit[tigHitIndAB].energy << ", array num: " << sortedEvt.tigHit[tigHitIndAB].core << ", address: " << sortedEvt.tigHit[tigHitIndAB]->GetAddress() << endl;

      numTigressABHits++;

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
        addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy);
        tigRate->Fill(sortedEvt.tigHit[tigHitIndAB].timeNs/pow(10,9));
        tigNum_time->Fill(sortedEvt.tigHit[tigHitIndAB].timeNs/pow(10,9),sortedEvt.tigHit[tigHitIndAB].core);
        //cout << "hit " << tigHitIndAB << ", Anum: " << sortedEvt.tigHit[tigHitIndAB].core << ", energy: " << sortedEvt.tigHit[tigHitIndAB].energy << endl;
        addE_ANum->Fill(sortedEvt.tigHit[tigHitIndAB].core, sortedEvt.tigHit[tigHitIndAB].energy);
        double theta = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180./PI;
        addE_theta->Fill(theta, sortedEvt.tigHit[tigHitIndAB].energy);
        //fill ring spectra
        addE_ring[getTIGRESSRing(theta)]->Fill(sortedEvt.tigHit[tigHitIndAB].energy);

        //check time-correlated TIP events
        uint16_t tipRingHP = 0;
        for (int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
          if((sortedEvt.csiHit[tipHitInd].detNum > 0)&&(sortedEvt.csiHit[tipHitInd].detNum <= NTIP)){
            numTipTigHits++;
            tipRingHP |= (1 << getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)); //update ring hitpattern
          }
        }

        for(int i=0; i<NTIPRING; i++){
          if(tipRingHP & (1 << i)){
            //there was a hit in TIP ring i
            tigE_tipRing->Fill(sortedEvt.tigHit[tigHitIndAB].energy,i);
          }
        }

        tigE_tipMult->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.header.numCsIHits);
      }
      
    }

    if((sortedEvt.header.numCsIHits>=0)||(sortedEvt.header.numTigABHits>=0)){
      tiptig_Amult->Fill(sortedEvt.header.numCsIHits,sortedEvt.header.numTigABHits);
    }

    //evaluate timing conditions
    //get fit times
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
      tipFitTimes[tipHitInd] = sortedEvt.csiHit[tipHitInd].timeNs;
    }

    //TIP-TIP timing, position, and energy
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
      for(int tipHitInd2 = tipHitInd+1; tipHitInd2 < sortedEvt.header.numCsIHits; tipHitInd2++){
        Double_t tDiff = tipFitTimes[tipHitInd] - tipFitTimes[tipHitInd2];
        tiptipFitT->Fill(tDiff);
        tipPos_tipPos->Fill(sortedEvt.csiHit[tipHitInd].detNum,sortedEvt.csiHit[tipHitInd2].detNum);
      }
    }
    
    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigABHits; tigHitIndAB++){

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){

        //TIP-TIG timing and energy
        for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
          Double_t tDiff = tipFitTimes[tipHitInd] - sortedEvt.tigHit[tigHitIndAB].timeNs;
          tipT_tigT_diff->Fill(tDiff);
          tigE_tipE->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.tigHit[tigHitIndAB].energy);
          tigE_tipTtigTdiff->Fill(tDiff,sortedEvt.tigHit[tigHitIndAB].energy);
        }

        //TIG-TIG timing, position, and energy
        for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigABHits; tigHitIndAB2++){
          if(sortedEvt.tigHit[tigHitIndAB2].energy > MIN_TIG_EAB){
            Double_t tDiff = sortedEvt.tigHit[tigHitIndAB].timeNs - sortedEvt.tigHit[tigHitIndAB2].timeNs;
            addT_addT->Fill(tDiff);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
            tigPos_tigPos->Fill(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB2].core);
          }
        }
      }
    }

    //count the number of protons or alphas
    for(int tipHitInd=0;tipHitInd<sortedEvt.header.numCsIHits;tipHitInd++){
      switch(getParticleTypePID(sortedEvt.csiHit[tipHitInd].PID,sortedEvt.csiHit[tipHitInd].energy, sortedEvt.csiHit[tipHitInd].detNum, gates)){ //see common.cxx
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
    
    if(evtNumProtons<=MAX_NUM_PARTICLE){
      if(evtNumAlphas<=MAX_NUM_PARTICLE){
        if(evtNumProtons+evtNumAlphas<=MAX_NUM_PARTICLE){

          for(int tigHitIndAB=0;tigHitIndAB<sortedEvt.header.numTigABHits;tigHitIndAB++){
            //cout << "energy: " << sortedEvt.tigHit[tigHitIndAB].energy << ", array num: " << sortedEvt.tigHit[tigHitIndAB].core << ", address: " << sortedEvt.tigHit[tigHitIndAB]->GetAddress() << endl;
            if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
              //TIGRESS PID separated addback energy
              Double_t angle = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180.0/PI;
              addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,0); //sum spectrum
              addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,getTIGRESSSegmentRing(angle)+1);
              double eDopp = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
              addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,0); //sum spectrum
              addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,getTIGRESSSegmentRing(angle)+1);
              for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigABHits; tigHitIndAB2++){
                if(sortedEvt.tigHit[tigHitIndAB2].energy > MIN_TIG_EAB){
                  addEaddE_xayp[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
                  addEaddE_xayp[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
                  double eDopp2 = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB2],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
                  addDoppaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(eDopp,eDopp2);
                  addDoppaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(eDopp2,eDopp); //symmetrized
                  addEaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,eDopp2);
                  addEaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,eDopp); //symmetrized
                }
              }
            }
          }

        }/*else{
          cout << "Entry " << jentry << " has too many charged particles (" << evtNumProtons+evtNumAlphas << ")!" << endl;
        }*/
      }/*else{
        cout << "Entry " << jentry << " has too many alphas (" << evtNumAlphas << ")!" << endl;
      }*/
    }/*else{
      cout << "Entry " << jentry << " has too many protons (" << evtNumProtons << ")!" << endl;
    }*/
    

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigressABHits << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP + TIGRESS hits: " << numTipTigHits << endl << endl;
  cout << endl << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  myfile->cd();

  TDirectory *tipdir = myfile->mkdir("TIP");
  tipdir->cd();
  tipList->Write();
  myfile->cd();

  TDirectory *tipPIDdir = myfile->mkdir("TIP PID Plots");
  tipPIDdir->cd();
  tipPIDList->Write();
  myfile->cd();

  TDirectory *tipPIDGatedir = myfile->mkdir("TIP PID Gates");
  tipPIDGatedir->cd();
  tipPIDGateList->Write();
  myfile->cd();

  TDirectory *timingdir = myfile->mkdir("Timing");
  timingdir->cd();
  timingList->Write();
  myfile->cd();

  TDirectory *tiptipdir = myfile->mkdir("TIP-TIP");
  tiptipdir->cd();
  tiptipList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  TDirectory *tiptigdir = myfile->mkdir("TIP-TIGRESS");
  tiptigdir->cd();
  tiptigList->Write();
  myfile->cd();

  TDirectory *tigPIDsepdir = myfile->mkdir("TIGRESS PID+Timing Separated");
  tigPIDsepdir->cd();
  tigPIDSepList->Write();
  myfile->cd();

  TDirectory *tigtigPIDsepdir = myfile->mkdir("TIGRESS-TIGRESS PID+Timing Separated");
  tigtigPIDsepdir->cd();
  tigtigPIDSepList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
  fclose(inp);
}
int main(int argc, char **argv)
{

  SortDiagnostics *mysort = new SortDiagnostics();

  char const *sfile;
  char const *outfile;
  printf("Starting SortDiagnosticsSMOL\n");

  // Input-chain-file, output-histogram-file
  if (argc == 1){
    cout << "Code sorts a bunch of diagnostic histograms for online TIP+TIGRESS data" << endl;
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
