//use the Makefile!

#define SortDiagnosticsS_cxx
#include "common.cxx"
#include "SortDiagnosticsSMOL.h"

using namespace std;

Int_t numTipRingHits[NTIPRING];
sorted_evt sortedEvt;
sorted_evt pastEvts[PAST_EVT_BUFSIZE];

void SortDiagnosticsS::SortData(const char *sfile, const char *outfile)
{
  Initialise();

  FILE *inp = fopen(sfile, "rb");
  printf("File %s opened\n", sfile);
  
  uint64_t sentries = 0U;
  fread(&sentries,sizeof(uint64_t),1,inp);
  
  uint8_t pastEvtPos = 0;
  uint8_t evtMixAllowed = 0;
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
    if(readSMOLEvent(inp,&sortedEvt)==0){
      cout << "ERROR: bad event data in entry " << jentry << "." << endl;
      exit(-1);
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

        tipRate->Fill(csiHitTime(&sortedEvt,tipHitInd)/pow(10,9));
        
        numTipRingHits[getTIPRing(sortedEvt.csiHit[tipHitInd].detNum)]++;
        numTipHits++;
        
        tip_E->Fill(sortedEvt.csiHit[tipHitInd].energy);
        tipEEvt += sortedEvt.csiHit[tipHitInd].energy;
        tip_pos->Fill(sortedEvt.csiHit[tipHitInd].detNum);
        tip_ring->Fill(getTIPRing(sortedEvt.csiHit[tipHitInd].detNum));
        tip_E_pos->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.csiHit[tipHitInd].detNum);
        tip_fittype->Fill(csiFitType(&sortedEvt.csiHit[tipHitInd]));

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

    for (int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
        tigE->Fill(sortedEvt.noABHit[tigHitInd].energy);
      }
    }

    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){

      //cout << "energy: " << sortedEvt.tigHit[tigHitIndAB].energy << ", array num: " << sortedEvt.tigHit[tigHitIndAB].core << ", address: " << sortedEvt.tigHit[tigHitIndAB]->GetAddress() << endl;

      numTigressABHits++;

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
        addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy);
        double eDopp = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
        addDopp->Fill(eDopp);
        /*//The following kills performance, uncomment only if needed
        for(Double_t beta=0.03;beta<0.05;beta+=0.0001){
          addDopp_dopp->Fill(beta,getEDoppFusEvapDirectBeta(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,beta,sortedEvt.csiHit,gates));
        }*/
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

    if((sortedEvt.header.numCsIHits>=0)||(sortedEvt.header.numTigHits>=0)){
      tiptig_Amult->Fill(sortedEvt.header.numCsIHits,sortedEvt.header.numTigHits);
    }

    //evaluate timing conditions
    //get fit times
    Double_t avgTIPHitT = 0;
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
      tipFitTimes[tipHitInd] = csiHitTime(&sortedEvt,tipHitInd);
      avgTIPHitT += csiHitTime(&sortedEvt,tipHitInd);
    }
    if(sortedEvt.header.numCsIHits > 0){
      avgTIPHitT /= sortedEvt.header.numCsIHits;
    }

    //TIP-TIP timing, position, and energy
    for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
      for(int tipHitInd2 = tipHitInd+1; tipHitInd2 < sortedEvt.header.numCsIHits; tipHitInd2++){
        Double_t tDiff = tipFitTimes[tipHitInd] - tipFitTimes[tipHitInd2];
        tiptipFitT->Fill(tDiff);
        tipPos_tipPos->Fill(sortedEvt.csiHit[tipHitInd].detNum,sortedEvt.csiHit[tipHitInd2].detNum);
        tipPos_tipPos->Fill(sortedEvt.csiHit[tipHitInd2].detNum,sortedEvt.csiHit[tipHitInd].detNum);
        tipRing_tipRing->Fill(getTIPRing(sortedEvt.csiHit[tipHitInd].detNum),getTIPRing(sortedEvt.csiHit[tipHitInd2].detNum));
        tipRing_tipRing->Fill(getTIPRing(sortedEvt.csiHit[tipHitInd2].detNum),getTIPRing(sortedEvt.csiHit[tipHitInd].detNum));
      }
    }
    
    for(int tigHitIndAB = 0; tigHitIndAB < sortedEvt.header.numTigHits; tigHitIndAB++){

      if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){

        //TIP-TIG timing and energy
        for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
          Double_t tDiff = tigHitTime(&sortedEvt,tigHitIndAB) - tipFitTimes[tipHitInd];
          tipT_tigT_diff->Fill(tDiff);
          tipTtigT_addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy,tigHitTime(&sortedEvt,tigHitIndAB)-avgTIPHitT);
          double eDopp = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
          Double_t ring = getTIGRESSRing(getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180.0/PI);
          if((ring!=2)&&(ring!=3)){
            tipTtigT_EDopp_no90->Fill(eDopp,tigHitTime(&sortedEvt,tigHitIndAB)-avgTIPHitT);
          }
          tipTtigT_EDopp->Fill(eDopp,tigHitTime(&sortedEvt,tigHitIndAB)-avgTIPHitT);
          tigE_tipE->Fill(sortedEvt.csiHit[tipHitInd].energy,sortedEvt.tigHit[tigHitIndAB].energy);
          tigE_tipTtigTdiff->Fill(tDiff,sortedEvt.tigHit[tigHitIndAB].energy);
        }

        //TIG-TIG timing, position, and energy
        for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
          if(sortedEvt.tigHit[tigHitIndAB2].energy > MIN_TIG_EAB){
            Double_t tDiff = tigHitTime(&sortedEvt,tigHitIndAB) - tigHitTime(&sortedEvt,tigHitIndAB2);
            addT_addT->Fill(tDiff);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy,sortedEvt.tigHit[tigHitIndAB2].energy);
            addE_addE->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
            tigPos_tigPos->Fill(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB2].core);
            tigTtigT_addE->Fill(sortedEvt.tigHit[tigHitIndAB].energy,fabs(tigHitTime(&sortedEvt,tigHitIndAB)-tigHitTime(&sortedEvt,tigHitIndAB2)));
            tigTtigT_addE->Fill(sortedEvt.tigHit[tigHitIndAB2].energy,fabs(tigHitTime(&sortedEvt,tigHitIndAB)-tigHitTime(&sortedEvt,tigHitIndAB2)));
          }
        }
      }
    }

    for(int tigHitInd = 0; tigHitInd < sortedEvt.header.numNoABHits; tigHitInd++){
      if(sortedEvt.noABHit[tigHitInd].energy > MIN_TIG_EAB){
        //TIP-TIG timing and energy (no addback)
        for(int tipHitInd = 0; tipHitInd < sortedEvt.header.numCsIHits; tipHitInd++){
          Double_t tDiffNoAB = noABHitTime(&sortedEvt,tigHitInd) - tipFitTimes[tipHitInd];
          tigE_tipTtigTdiffNoAB->Fill(tDiffNoAB,sortedEvt.tigHit[tigHitInd].energy);
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

          for(int tigHitIndAB=0;tigHitIndAB<sortedEvt.header.numTigHits;tigHitIndAB++){
            //cout << "energy: " << sortedEvt.tigHit[tigHitIndAB].energy << ", array num: " << sortedEvt.tigHit[tigHitIndAB].core << ", address: " << sortedEvt.tigHit[tigHitIndAB]->GetAddress() << endl;
            if(sortedEvt.tigHit[tigHitIndAB].energy > MIN_TIG_EAB){
              //TIGRESS PID separated addback energy
              Double_t angle = getTigVector(sortedEvt.tigHit[tigHitIndAB].core,sortedEvt.tigHit[tigHitIndAB].seg).Theta()*180.0/PI;
              addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,0); //sum spectrum
              addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,getTIGRESSSegmentRing(angle)+1);
              double eDopp = getEDoppFusEvapDirect(&sortedEvt.tigHit[tigHitIndAB],sortedEvt.header.numCsIHits,sortedEvt.csiHit,gates);
              addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,0); //sum spectrum
              addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,getTIGRESSSegmentRing(angle)+1);
              for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < sortedEvt.header.numTigHits; tigHitIndAB2++){
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
              if(evtMixAllowed == 2){
                uint8_t evtMixEvt = pastEvtPos+5;
                if(evtMixEvt >= PAST_EVT_BUFSIZE){
                  evtMixEvt = 0;
                }
                for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < pastEvts[evtMixEvt].header.numTigHits; tigHitIndAB2++){
                  addEaddE_xayp_evtmix[evtNumProtons][evtNumAlphas]->Fill(sortedEvt.tigHit[tigHitIndAB].energy,pastEvts[evtMixEvt].tigHit[tigHitIndAB2].energy);
                  addEaddE_xayp_evtmix[evtNumProtons][evtNumAlphas]->Fill(pastEvts[evtMixEvt].tigHit[tigHitIndAB2].energy,sortedEvt.tigHit[tigHitIndAB].energy); //symmetrized
                }
              }
            }
          }

          memcpy(&pastEvts[pastEvtPos],&sortedEvt,sizeof(sortedEvt));
          pastEvtPos++;
          if(pastEvtPos >= PAST_EVT_BUFSIZE){
            if(evtMixAllowed < 2){
              evtMixAllowed++;
            }
            pastEvtPos = 0;
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

  TDirectory *tipPIDdir = myfile->mkdir("TIP_PID_Plots");
  tipPIDdir->cd();
  tipPIDList->Write();
  myfile->cd();

  TDirectory *tipPIDGatedir = myfile->mkdir("TIP_PID_Gates");
  tipPIDGatedir->cd();
  tipPIDGateList->Write();
  myfile->cd();

  TDirectory *timingdir = myfile->mkdir("Timing");
  timingdir->cd();
  timingList->Write();
  myfile->cd();

  TDirectory *tiptipdir = myfile->mkdir("TIP_TIP");
  tiptipdir->cd();
  tiptipList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS_TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  TDirectory *tiptigdir = myfile->mkdir("TIP_TIGRESS");
  tiptigdir->cd();
  tiptigList->Write();
  myfile->cd();

  TDirectory *tigPIDsepdir = myfile->mkdir("TIGRESS_PID_and_Timing_Separated");
  tigPIDsepdir->cd();
  tigPIDSepList->Write();
  myfile->cd();

  TDirectory *tigtigPIDsepdir = myfile->mkdir("TIGRESS_TIGRESS_PID_and_Timing_Separated");
  tigtigPIDsepdir->cd();
  tigtigPIDSepList->Write();
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
