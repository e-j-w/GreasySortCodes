//g++ SortDiagnostics.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData

#define SortDiagnostics_cxx
#include "common.h"
#include "SortDiagnostics.h"

using namespace std;

float lastTIPHitT[NTIP]; //stores the last hit time for each detector
Int_t numTipRingPileupHits[NTIPRING], numTipRingHits[NTIPRING];

bool suppTig = false;
bool suppAdd = false;

void SortDiagnostics::SortData(char const *afile, char const *calfile, char const *outfile)
{
  Initialise();

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()){
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if(AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }else{
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  TTip *tip = 0;
  if (AnalysisTree->FindBranch("TTip")){
    AnalysisTree->SetBranchAddress("TTip", &tip);
  }else{
    cout << "Branch 'TTip' not found! TTip variable is NULL pointer" << endl;
  }
  
  unsigned long int numTigressABHits = 0;
  unsigned long int numTipHits = 0;
  unsigned long int numTipPileupHits = 0;
  unsigned long int numTipTigHits = 0;

  unsigned long int numTigressEvt = 0;
  unsigned long int numTigPileupHits = 0;
  unsigned long int numTipEvt = 0;
  unsigned long int numTipTigEvt = 0;

  Double_t tipFitTimes[MAXNUMTIPHIT];

  //Defining Pointers
  TTigressHit *tig_hit, *add_hit, *add_hit2;
  TTipHit *tip_hit, *tip_hit2;
  const std::vector<Short_t> *wf; //for CsI waveform

  double_t tTipCFD;
  double tipEEvt;
  Int_t evtNumProtons, evtNumAlphas;
  Int_t evtNumProtonsDetSumGate, evtNumAlphasDetSumGate;
  Int_t evtNumHorDetSumGate;
  memset(lastTIPHitT,0,sizeof(lastTIPHitT));

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  tigress->SetSuppressionCriterion(ExptSuppression);

  /*cout << "TIGRESS positions: " << endl;
  for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
      TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
      cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
  }*/
  

  printf("\nSorting analysis events...\n");
  for(Long64_t jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
      continue;
    }

    if(tip && tip->GetMultiplicity()>MAXNUMTIPHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIP hits (" << tip->GetMultiplicity() << ")!" << endl;
      continue;
    }
    if(tigress && tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
      cout << "Ignoring entry " << jentry << " as it has too many TIGRESS hits (" << tigress->GetAddbackMultiplicity() << ")!" << endl;
      continue;
    }
    
    evtNumProtons=0;
    evtNumAlphas=0;
    evtNumProtonsDetSumGate=0;
    evtNumAlphasDetSumGate=0;
    evtNumHorDetSumGate=0;
    
    if(tip){
      
      numTipEvt++;
      if(tip->GetMultiplicity()>=0){
        tip_mult->Fill(tip->GetMultiplicity());
      }

      tipEEvt = 0.;
      for (int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){

        tip_hit = tip->GetTipHit(tipHitInd);
        if((tip_hit->GetTipChannel() > 0)&&(tip_hit->GetTipChannel() <= NTIP)){

          tipRate->Fill(tip_hit->GetTime()/pow(10,9));
          
          //cout << "TIP hit " << tipHitInd << ", K: " << tip_hit->GetKValue() << endl;
          if(tip_hit->GetKValue() != noPileupKValue){
            //pileup
            numTipRingPileupHits[getTIPRing(tip_hit->GetTipChannel())]++;
            numTipPileupHits++;
            continue;
          }
          
          numTipRingHits[getTIPRing(tip_hit->GetTipChannel())]++;
          numTipHits++;

          if(lastTIPHitT[tip_hit->GetTipChannel()-1]==0.){
            lastTIPHitT[tip_hit->GetTipChannel()-1] = tip_hit->GetTime();
          }else{
            //cout << "val: " << tip_hit->GetTime() - lastTIPHitT[tip_hit->GetTipChannel()-1] << endl;
            lastTIPHitT[tip_hit->GetTipChannel()-1] = tip_hit->GetTime();
          }

          wf = tip_hit->GetWaveform();
          tip_wfrmsize->Fill(wf->size());
          tip_fittype->Fill(tip_hit->GetFitType());
          tTipCFD = tip_hit->GetTime() - (((tip_hit->GetTimeStamp()) + gRandom->Uniform()) * tip_hit->GetTimeStampUnit());
          tTipCFD += 1000.; //offset by pretrigger
          //cout << "CFD time: " << tTipCFD << ", fit time: " << tTipFit << endl;
          
          tip_E->Fill(tip_hit->GetEnergy());
          tipEEvt += tip_hit->GetEnergy();
          tip_pos->Fill(tip_hit->GetTipChannel());
          tip_ring->Fill(getTIPRing(tip_hit->GetTipChannel()));
          tip_CFDFitDiff->Fill(getTipFitTime(tip_hit,tip_waveform_pretrigger) - tTipCFD);
          tip_E_pos->Fill(tip_hit->GetEnergy(),tip_hit->GetTipChannel());

          //PID related stuff
          if(tip_hit->GetPID() >= 0){ //PID was found
            tip_E_PID_Sum->Fill(tip_hit->GetEnergy(),tip_hit->GetPID());
            //cout << "Energy: " << tip_hit->GetEnergy() << ", PID: " << tip_hit->GetPID() << endl;
            tip_E_PID[tip_hit->GetTipChannel()-1]->Fill(tip_hit->GetEnergy(),tip_hit->GetPID());
            tip_E_PID_Ring[getTIPRing(tip_hit->GetTipChannel())]->Fill(tip_hit->GetEnergy(),tip_hit->GetPID());
          }
        }
        
      }
      tip_Etot->Fill(tipEEvt);
    }

    if(tigress){
      numTigressEvt++;

      for (int tigHitInd = 0; tigHitInd < tigress->GetMultiplicity(); tigHitInd++){
        tig_hit = tigress->GetTigressHit(tigHitInd);
        suppTig = tig_hit->BGOFired();

        //cout << "TIGRESS hit " << tigHitInd << ", K: " << tig_hit->GetKValue() << endl;
        if (noPileupKValue != tig_hit->GetKValue()){
          numTigPileupHits++;
        }

        //cout << "TIGRESS hit " << tigHitInd << ", arrayNum: " << tig_hit->GetArrayNumber() << ", charge: " << tig_hit->GetCharge() << endl;
        if(!suppTig && tig_hit->GetEnergy() > MIN_TIG_EAB){
          tigE->Fill(tig_hit->GetEnergy());
          tigE_ANum->Fill(tig_hit->GetArrayNumber(), tig_hit->GetEnergy());
          tigChan->Fill(tig_hit->GetArrayNumber());
          tigRate->Fill(tig_hit->GetTime()/pow(10,9));
          for(int bgoInd=0; bgoInd < tigress->GetBGOMultiplicity(); bgoInd++){
            if((tigress->GetBGO(bgoInd).GetDetector() == tig_hit->GetDetector()) && (tigress->GetBGO(bgoInd).GetEnergy() > 0.)){
              tigT_bgoT_supp->Fill(tig_hit->GetCfd() - tigress->GetBGO(bgoInd).GetCfd());
            }
          }
        }
        if(tig_hit->GetEnergy() > MIN_TIG_EAB){
          tigE_unsupp->Fill(tig_hit->GetEnergy());

          bgo_mult->Fill(tigress->GetBGOMultiplicity());
          for(int bgoInd=0; bgoInd < tigress->GetBGOMultiplicity(); bgoInd++){
            bgo_det->Fill(tigress->GetBGO(bgoInd).GetDetector());
            if((tigress->GetBGO(bgoInd).GetDetector() == tig_hit->GetDetector()) && (tigress->GetBGO(bgoInd).GetEnergy() > 0.)){
              tigT_bgoT->Fill(tig_hit->GetCfd() - tigress->GetBGO(bgoInd).GetCfd());
            }
          }
        }

      }

      for (int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
        add_hit = tigress->GetAddbackHit(tigHitIndAB);
        //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
        if (!add_hit->BGOFired() && add_hit->GetEnergy() > MIN_TIG_EAB){

          numTigressABHits++;
	  

          addE->Fill(add_hit->GetEnergy());
          //cout << "hit " << tigHitIndAB << ", Anum: " << add_hit->GetArrayNumber() << ", energy: " << add_hit->GetEnergy() << endl;
          addE_ANum->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
          float theta = add_hit->GetPosition().Theta()*180./PI;
          addE_theta->Fill(theta, add_hit->GetEnergy());
          //fill ring spectra
          addE_ring[getTIGRESSRing(theta)]->Fill(add_hit->GetEnergy());

          //check time-correlated TIP events
          if(tip){
            
            uint16_t tipRingHP = 0;
            
            for (int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){

              tip_hit = tip->GetTipHit(tipHitInd);
              double tTipFit = getTipFitTime(tip_hit,tip_waveform_pretrigger);
              if((tip_hit->GetTipChannel() > 0)&&(tip_hit->GetTipChannel() <= NTIP)){
                numTipTigHits++;
                tipRingHP |= (1 << getTIPRing(tip_hit->GetTipChannel())); //update ring hitpattern
              }
            }

            for(int i=0; i<NTIPRING; i++){
              if(tipRingHP & (1 << i)){
                //there was a hit in TIP ring i
                tigE_tipRing->Fill(add_hit->GetEnergy(),i);
              }
            }

            tigE_tipMult->Fill(add_hit->GetEnergy(),tip->GetMultiplicity());

          }
          
        }else if(add_hit->BGOFired() && add_hit->GetEnergy() > MIN_TIG_EAB){
          addE_rej->Fill(add_hit->GetEnergy());
        }
      }
    }

    if(tip && tigress){
      numTipTigEvt++;

      //check whether event passes all timing gates (TIG-TIG, TIP-TIP, TIP-TIG)
      //also rejects pileup
      uint64_t passedtimeGate = passesTimeGate(tigress,tip,1,2);

      //ensure full TIP-TIGRESS timing condition is met
      //ie. do not accept hits with partial timing condition match
      if(!(passedtimeGate&(1ULL<<TIPTIGFLAG))){
        passedtimeGate = 0;
      }

      //hit multiplicities
      if((tip->GetMultiplicity()>=0)||(tigress->GetAddbackMultiplicity()>=0)){
        Int_t tigSuppMult=0;
        Int_t tigSuppTimingMult=0;
        Int_t tipTimingMult=0;
        for (int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
          add_hit = tigress->GetAddbackHit(tigHitIndAB);
          suppAdd = add_hit->BGOFired();
          if(!suppAdd && add_hit->GetEnergy() > MIN_TIG_EAB){
            tigSuppMult++;
            if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
              tigSuppTimingMult++;
            }
          }
        }
        for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
           if(passedtimeGate&(1ULL<<tipHitInd)){
             tipTimingMult++;
           }
        }
        tiptig_multSupp->Fill(tip->GetMultiplicity(),tigSuppMult);
        if((tipTimingMult>0)||(tigSuppTimingMult>0)){
          tiptig_multSuppTiming->Fill(tipTimingMult,tigSuppTimingMult);
        }
      }
      if((tip->GetMultiplicity()>=0)||(tigress->GetMultiplicity()>=0)){
        tiptig_mult->Fill(tip->GetMultiplicity(),tigress->GetMultiplicity());
      }

      //evaluate timing conditions
      //get fit times
      for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
        tip_hit = tip->GetTipHit(tipHitInd);
        tipFitTimes[tipHitInd] = getTipFitTime(tip_hit,tip_waveform_pretrigger);
      }

      //TIP-TIP timing
      for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
        numTipHits++;
        tip_hit = tip->GetTipHit(tipHitInd);
        for(int tipHitInd2 = tipHitInd+1; tipHitInd2 < tip->GetMultiplicity(); tipHitInd2++){
          tip_hit2 = tip->GetTipHit(tipHitInd2);
          Double_t tDiff = tipFitTimes[tipHitInd] - tipFitTimes[tipHitInd2];
          tiptipFitT->Fill(tDiff);
          //tiptipFitT->Fill(tip_hit->GetTime() - tip_hit2->GetTime());
          if(passedtimeGate&(1ULL<<tipHitInd)){
            if(passedtimeGate&(1ULL<<tipHitInd2)){
              tiptipFitTPassed->Fill(tDiff);
            }
          }
        }
      }
      
      for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
        add_hit = tigress->GetAddbackHit(tigHitIndAB);
        suppAdd = add_hit->BGOFired();
        //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
        if (!suppAdd && add_hit->GetEnergy() > MIN_TIG_EAB){

          //TIP-TIG timing and energy
          for(int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){
            tip_hit = tip->GetTipHit(tipHitInd);
            Double_t tDiff = tipFitTimes[tipHitInd] - add_hit->GetTime();
            Double_t tDiffCFD = tip_hit->GetTime() - add_hit->GetTime();
            tipT_tigT_diff->Fill(tDiff);
            tipTCFD_tigT_diff->Fill(tDiffCFD);
            if(passedtimeGate&(1ULL<<tipHitInd)){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                tipT_tigT_diffPassed->Fill(tDiff);
                tipTCFD_tigT_diffPassed->Fill(tDiffCFD);
              }
            }
            tigE_tipTtigTdiff->Fill(tDiff,add_hit->GetEnergy());
          }

          //TIG-TIG timing and energy
          for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < tigress->GetAddbackMultiplicity(); tigHitIndAB2++){
            add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
            suppAdd = add_hit2->BGOFired();
            if(!suppAdd && add_hit2->GetEnergy() > MIN_TIG_EAB){
              Double_t tDiff = add_hit->GetTime() - add_hit2->GetTime();
              addT_addT->Fill(tDiff);
              addE_addE->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
              addE_addE->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                  addT_addTPassed->Fill(tDiff);
                  addE_addE_tgPassed->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                  addE_addE_tgPassed->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
                }
              }
            }
          }
        }
      }

      //count the number of protons or alphas
      for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
        if(passedtimeGate&(1ULL<<tipHitInd)){
          tip_hit = tip->GetTipHit(tipHitInd);
          numTipHits++;

          switch(getParticleType(tip_hit,gates)){ //see common.cxx
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
      }
      
      if(evtNumProtons<=MAX_NUM_PARTICLE){
        if(evtNumAlphas<=MAX_NUM_PARTICLE){
          if(evtNumProtons+evtNumAlphas<=MAX_NUM_PARTICLE){

            for(int tigHitIndAB=0;tigHitIndAB<tigress->GetAddbackMultiplicity();tigHitIndAB++){
              if(passedtimeGate&(1ULL<<(tigHitIndAB+MAXNUMTIPHIT))){
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                suppAdd = add_hit->BGOFired();
                //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
                if(!suppAdd && add_hit->GetEnergy() > MIN_TIG_EAB){
                  //TIGRESS PID separated addback energy
                  double_t thetaDeg = add_hit->GetPosition().Theta()*180./PI;
                  addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(add_hit->GetEnergy(),0); //sum spectrum
                  addE_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(add_hit->GetEnergy(),getTIGRESSRing(thetaDeg)+1);
                  double eDopp = getEDoppFusEvap(add_hit,tip,passedtimeGate,gates);
                  addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,0); //sum spectrum
                  addDopp_xayp_ring[evtNumProtons][evtNumAlphas]->Fill(eDopp,getTIGRESSRing(thetaDeg)+1);
                  for(int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < tigress->GetAddbackMultiplicity(); tigHitIndAB2++){
                    if(passedtimeGate&(1ULL<<(tigHitIndAB2+MAXNUMTIPHIT))){
                      add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
                      suppAdd = add_hit2->BGOFired();
                      if(!suppAdd && add_hit2->GetEnergy() > MIN_TIG_EAB){
                        double eDopp2 = getEDoppFusEvap(add_hit2,tip,passedtimeGate,gates);
                        addEaddE_xayp[evtNumProtons][evtNumAlphas]->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                        addEaddE_xayp[evtNumProtons][evtNumAlphas]->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
                        addDoppaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(eDopp,eDopp2);
                        addDoppaddDopp_xayp[evtNumProtons][evtNumAlphas]->Fill(eDopp2,eDopp); //symmetrized
                      }
                    }
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

    }
    

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigressABHits << endl;
  cout << "Number of TIGRESS pileup hits: " << numTigPileupHits << " (" << (float)(numTigPileupHits*100)/float(numTigPileupHits + numTigressABHits) << " %)" << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP pileup hits: " << numTipPileupHits << " (" << (float)(numTipPileupHits*100)/float(numTipPileupHits + numTipHits) << " %)" << endl;
  cout << "Number of TIP + TIGRESS hits: " << numTipTigHits << endl << endl;
  /*cout << "Number of TIGRESS events: " << numTigressEvt << endl;
  cout << "Number of TIP events: " << numTipEvt << endl;
  cout << "Number of TIP + TIGRESS events: " << numTipTigEvt << endl << endl;*/
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  myfile->cd();

  TDirectory *tigbgodir = myfile->mkdir("TIGRESS-BGO");
  tigbgodir->cd();
  tigbgoList->Write();
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
  analysisfile->Close();
}
int main(int argc, char **argv)
{

  SortDiagnostics *mysort = new SortDiagnostics();

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting SortDiagnostics\n");

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
    cout << "Arguments: SortDiagnostics analysis_tree calibration_file output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }else if(argc == 2){
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile); 
  }else if(argc == 3){
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else if(argc == 4){
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }else{
    printf("ERROR: too many arguments!\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
