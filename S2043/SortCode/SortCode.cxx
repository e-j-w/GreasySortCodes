//g++ SortCode.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData

#define Sortcode_cxx
#include "SortCode.h"

using namespace std;

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

//function defines the Compton suppression timing window for TIGRESS
bool S2043Suppression(TDetectorHit* tig, TBgoHit& bgo)
{
   Int_t dCfd = static_cast<TTigressHit*>(tig)->GetCfd() - bgo.GetCfd();
   return ((dCfd > -50 && dCfd < 350) && (tig->GetDetector() == bgo.GetDetector()) && (bgo.GetEnergy() > 0));
}


bool gate1D(Double_t value, Double_t min, Double_t max)
{
  if (min < value && value < max)
    return true;
  else
    return false;
}

Int_t getTIPRing(const Int_t tipPosition){
  if((tipPosition < 1)||(tipPosition > 128)){
    cout << "Invalid TIP position: " << tipPosition << endl;
    return 0;
  }else if(tipPosition <= 4){
    return 0;
  }else if(tipPosition <= 10){
    return 1;
  }else if(tipPosition <= 22){
    return 2;
  }else if(tipPosition <= 38){
    return 3;
  }else if(tipPosition <= 58){
    return 4;
  }else if(tipPosition <= 76){
    return 5;
  }else if(tipPosition <= 94){
    return 6;
  }else if(tipPosition <= 108){
    return 7;
  }else if(tipPosition <= 120){
    return 8;
  }else{
    return 9;
  }
}

Int_t getTIGRESSRing(const float theta){
  if((theta < 0.0f)||(theta > 180.0f)){
    cout << "Invalid TIGRESS angle: " << theta << endl;
    return 0;
  }else if(theta <= 45.0f){
    return 0;
  }else if(theta <= 66.0f){
    return 1;
  }else if(theta <= 90.0f){
    return 2;
  }else if(theta <= 112.0f){
    return 3;
  }else if(theta <= 135.0f){
    return 4;
  }else{
    return 5;
  }
}


double pi = TMath::Pi();
double tigtigT[2] = {-100, 100}; // TIGRESS - TIGRESS Timing. Change !!!
double tiptigT[2] = {-575, -396}; // TIP - TIGRESS Fit Timing. Change !!!
double tiptigCFDT[2] = {-3000, 3000}; // TIP - TIGRESS CFD Timing. Change !!!
Int_t tip_waveform_pretrigger = 250;
float lastTIPHitT[NTIP]; //stores the last hit time for each detector
Int_t numTipRingPileupHits[NTIPRING], numTipRingHits[NTIPRING];


bool suppTig = false;
bool suppAdd = false;

void SortCode::SortData(char const *afile, char const *calfile, char const *outfile)
{
  Initialise();

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen())
  {
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long int analentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress")){
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }
  else
  {
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
  unsigned long int numTipEvt = 0;
  unsigned long int numTipTigEvt = 0;

  //Defining Pointers
  TTigressHit *tig_hit, *add_hit, *add_hit2;
  TTipHit *tip_hit;
  const std::vector<Short_t> *wf; //for CsI waveform

  double_t tTipCFD, tTipFit;
  double_t tipPID;
  Int_t evtNumProtons, evtNumAlphas;
  Int_t evtNumProtonsDetSumGate, evtNumAlphasDetSumGate;
  Int_t evtNumHorDetSumGate;
  memset(lastTIPHitT,0,sizeof(lastTIPHitT));

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  tigress->SetSuppressionCriterion(S2043Suppression);

  /*cout << "TIGRESS positions: " << endl;
  for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
      TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
      cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
  }*/
  

  printf("\nSorting analysis events...\n");
  for (int jentry = 0; jentry < analentries; jentry++){

    if(AnalysisTree->GetEntry(jentry) == 0){
      //entry not read successfully
      cout << "WARNING: could not read entry " << jentry << endl; 
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
      if(tip->GetMultiplicity()>=2){
        for(int i=1;i<tip->GetMultiplicity();i++){
          tiptipT->Fill(tip->GetTipHit(i-1)->GetTime() - tip->GetTipHit(i)->GetTime());
        }
      }

      //adding a tip-tip matrix to troubleshoot the multiplicty 0 events that were appearing
      if (2 == tip->GetMultiplicity()){
	      tipTipDetectors->Fill(tip->GetTipHit(0)->GetTipChannel(),tip->GetTipHit(1)->GetTipChannel());
      }	

      for (int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){

        tip_hit = tip->GetTipHit(tipHitInd);
        if(tip_hit->GetKValue() != 700){
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
          tip_lastEvtTDiff_pos->Fill((tip_hit->GetTime() - lastTIPHitT[tip_hit->GetTipChannel()-1])/1000.,tip_hit->GetTipChannel()); //in us
          lastTIPHitT[tip_hit->GetTipChannel()-1] = tip_hit->GetTime();
        }

        wf = tip_hit->GetWaveform();
        tip_wfrmsize->Fill(wf->size());
        TPulseAnalyzer pulse;
        pulse.SetData(*wf, 0);
        if(wf->size() > 50){
          double baseline = 0;
          double max, amplitude;
          for (unsigned int i = 1; i < wf->size()-2; i++){
            if (i < 51)
              baseline += wf->at(i);
            if (i == 1)
              max = wf->at(i);
            if (max < wf->at(i))
              max = wf->at(i);
            //cout << "sample " << i << ": " << wf->at(i) << endl;
          }
          baseline = baseline / 50.0;
          amplitude = (abs(max - baseline));
          tip_E_waveformAmp->Fill(tip_hit->GetEnergy(),amplitude);
          //cout << "amplitude: " << amplitude << endl;
          //getc(stdin);
          tipPID = pulse.CsIPID();
          if(tipPID > -1000.0) //in (modified) GRSISort, failed fits given value of -1000.0
            tipPID += 100.;
          tTipFit = pulse.CsIt0() * 10.;
          
        }
        tTipCFD = tip_hit->GetTime() - (((tip_hit->GetTimeStamp()) + gRandom->Uniform()) * tip_hit->GetTimeStampUnit());
        tTipCFD += 1000.; //offset by pretrigger
        //cout << "CFD time: " << tTipCFD << ", fit time: " << tTipFit << endl;
        
        tip_E->Fill(tip_hit->GetEnergy());
        tip_pos->Fill(tip_hit->GetTipChannel());
        tip_ring->Fill(getTIPRing(tip_hit->GetTipChannel()));
        tip_CFDFitDiff->Fill(tTipFit - tTipCFD);
        if(tip_hit->GetEnergy()>100.){
          tip_CFDFitDiffHighE->Fill(tTipFit - tTipCFD);
          /*TCanvas *c1 = new TCanvas("c1","Histogram",200,10,1200,1000);
          //pulse.DrawT0fit();
          pulse.DrawCsIFit();
          theApp->Run(kTRUE);*/
        }
        tip_E_pos->Fill(tip_hit->GetEnergy(),tip_hit->GetTipChannel());

        //PID related stuff
        if(tipPID>=0){ //PID was found
          tip_E_PID_Sum->Fill(tip_hit->GetEnergy(),tipPID);

          if((tip_hit->GetTipChannel() > 0)&&(tip_hit->GetTipChannel() <= NTIP)){
            tip_E_PID[tip_hit->GetTipChannel()-1]->Fill(tip_hit->GetEnergy(),tipPID);

            //count the number of each particle type
            if(alphaPosCut[tip_hit->GetTipChannel()-1]->IsInside(tip_hit->GetEnergy(),tipPID)){ //check if the hit is inside the alpha gate for the given TIP detector
              evtNumAlphas++;
            }else if(protonPosCut[tip_hit->GetTipChannel()-1]->IsInside(tip_hit->GetEnergy(),tipPID)){ //check if the hit is inside the proton gate for the given TIP detector
              evtNumProtons++;
            }
            if(alphaPosCutSum->IsInside(tip_hit->GetEnergy(),tipPID)){ //check if the hit is inside the alpha gate (common gate for each detector)
              evtNumAlphasDetSumGate++;
            }else if(protonPosCutSum->IsInside(tip_hit->GetEnergy(),tipPID)){ //check if the hit is inside the proton gate (common gate for each detector)
              evtNumProtonsDetSumGate++;
            }else if(horPosCutSum->IsInside(tip_hit->GetEnergy(),tipPID)){ //check if hit is inside horizontal pid gate (common gate for each detector)
	      evtNumHorDetSumGate++;
	    }
          }
        }

        /*cout << "Drawing waveform" << endl;
        TCanvas *c1 = new TCanvas("c1","Waveform",200,10,1200,1000);
        TH1I *whist = pulse.GetWaveHist();
        char histName[256];
        whist->GetXaxis()->SetTitle("Sample Number");
        whist->SetStats(0);
        whist->Draw();
        theApp->Run(kTRUE);*/
        
      }
      if(tip->GetMultiplicity()>0){
        tip_particle_mult->Fill(evtNumProtons,evtNumAlphas);
      }
    }

    if(tigress){
      numTigressEvt++;

      if (2 == tigress->GetMultiplicity()){
	tigTigDetectors->Fill(tigress->GetTigressHit(0)->GetArrayNumber(),tigress->GetTigressHit(1)->GetArrayNumber());
      }

      if (2 == tigress->GetAddbackMultiplicity()){
	tigTigCloverDetectors->Fill(tigress->GetAddbackHit(0)->GetArrayNumber()/4,tigress->GetAddbackHit(1)->GetArrayNumber()/4);
      }


      for (int tigHitInd = 0; tigHitInd < tigress->GetMultiplicity(); tigHitInd++){
        tig_hit = tigress->GetTigressHit(tigHitInd);
        suppTig = tig_hit->BGOFired();
        if(!suppTig && tig_hit->GetEnergy() > 15){
          tigE->Fill(tig_hit->GetEnergy());
          tigE_ANum->Fill(tig_hit->GetArrayNumber(), tig_hit->GetEnergy());
          for(int bgoInd=0; bgoInd < tigress->GetBGOMultiplicity(); bgoInd++){
            if((tigress->GetBGO(bgoInd).GetDetector() == tig_hit->GetDetector()) && (tigress->GetBGO(bgoInd).GetEnergy() > 0.)){
              tigT_bgoT_supp->Fill(tig_hit->GetCfd() - tigress->GetBGO(bgoInd).GetCfd());
            }
          }
        }
        if(tig_hit->GetEnergy() > 15){
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
        suppAdd = add_hit->BGOFired();
        //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
        if (!suppAdd && add_hit->GetEnergy() > 15){

          numTigressABHits++;
          addE->Fill(add_hit->GetEnergy());
          //cout << "hit " << tigHitIndAB << ", Anum: " << add_hit->GetArrayNumber() << ", energy: " << add_hit->GetEnergy() << endl;
          addE_ANum->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
          float theta = add_hit->GetPosition().Theta()*180./PI;
          addE_theta->Fill(theta, add_hit->GetEnergy());
          num_addr->Fill(add_hit->GetArrayNumber(),add_hit->GetAddress());
          //fill ring spectra
          addE_ring[getTIGRESSRing(theta)]->Fill(add_hit->GetEnergy());

          //create gamma spectra for the two lampshades and the corona
          if (add_hit->GetArrayNumber() < 16){
          addE_ANum_Downstream->Fill(add_hit->GetEnergy());
          }
          else if (add_hit->GetArrayNumber() > 47){
          addE_ANum_Upstream->Fill(add_hit->GetEnergy());
          }
          else{
          addE_ANum_Corona->Fill(add_hit->GetEnergy());
          }

          //check time-correlated TIP events
          if(tip){

            uint16_t tipRingHP = 0;
            bool tipHitInCoinc = false; //whether or not there is any TIP hit in time coinc with this TIGRESS hit

            for (int tipHitInd = 0; tipHitInd < tip->GetMultiplicity(); tipHitInd++){

              tip_hit = tip->GetTipHit(tipHitInd);
              numTipTigHits++;
              tipRingHP |= (1 << getTIPRing(tip_hit->GetTipChannel())); //update ring hitpattern
              wf = tip_hit->GetWaveform();
              tip_wfrmsize->Fill(wf->size());
              TPulseAnalyzer pulse;
              pulse.SetData(*wf, 0);
              if(wf->size() > 50){
                tTipFit = pulse.CsIt0() * 10.;
                tTipFit += ((tip_hit->GetTimeStamp()) + gRandom->Uniform()) * tip_hit->GetTimeStampUnit();
                tTipFit -= tip_waveform_pretrigger * 10.;
                //cout << "TIP fit time: " << tTipFit << "TIGRESS time: " << add_hit->GetTime() << ", diff: " << tTipFit - add_hit->GetTime() << endl;
              }
              tipT_tigT_diff->Fill(tTipFit - add_hit->GetTime());
              tipTCFD_tigT_diff->Fill(tip_hit->GetTime() - add_hit->GetTime());

              //CFD time gated TIP
              if(gate1D(tip_hit->GetTime() - add_hit->GetTime(), tiptigCFDT[0], tiptigCFDT[1])){
                tipHitInCoinc = true;
                tipE_tigtg->Fill(tip_hit->GetEnergy());
                tipT_tigT_difftg->Fill(tTipFit - add_hit->GetTime());
                tipTCFD_tigT_difftg->Fill(tip_hit->GetTime() - add_hit->GetTime());
              }
              /*//fit time gated TIP
              if(gate1D(tTipFit - add_hit->GetTime(), tiptigT[0], tiptigT[1])){
                tipHitInCoinc = true;
                tipE_tigtg->Fill(tip_hit->GetEnergy());
                tipT_tigT_difftg->Fill(tTipFit - add_hit->GetTime());
                tipTCFD_tigT_difftg->Fill(tip_hit->GetTime() - add_hit->GetTime());
              }*/
              
            }

            //fit time gated TIGRESS
            if(tipHitInCoinc){
              //don't double count TIGRESS hits
              tigE_TIPtg->Fill(add_hit->GetEnergy());

              //particle gated TIGRESS spectra
              if(evtNumAlphas==2 && evtNumProtons==0){
                tigE_0p2aGate->Fill(add_hit->GetEnergy());
              }else if(evtNumAlphas==0 && evtNumProtons==2){
                tigE_2p0aGate->Fill(add_hit->GetEnergy());
              }
              if(evtNumAlphasDetSumGate==2 && evtNumProtonsDetSumGate==0){
                tigE_0p2aDetSumGate->Fill(add_hit->GetEnergy());

		//separate the spectra based on array number
		if (add_hit->GetArrayNumber() < 16){
			tigE_0p2aDetSumGate_Downstream->Fill(add_hit->GetEnergy());		
		}else if (add_hit->GetArrayNumber() > 47){
			tigE_0p2aDetSumGate_Upstream->Fill(add_hit->GetEnergy());		
		}else{
			tigE_0p2aDetSumGate_Corona->Fill(add_hit->GetEnergy());		
		}

              }else if(evtNumAlphasDetSumGate==0 && evtNumProtonsDetSumGate==2){
                tigE_2p0aDetSumGate->Fill(add_hit->GetEnergy());
		
		//separate the spectra based on array number
		if (add_hit->GetArrayNumber() < 16){
			tigE_2p0aDetSumGate_Downstream->Fill(add_hit->GetEnergy());		
		}else if (add_hit->GetArrayNumber() > 47){
			tigE_2p0aDetSumGate_Upstream->Fill(add_hit->GetEnergy());		
		}else{
			tigE_2p0aDetSumGate_Corona->Fill(add_hit->GetEnergy());		
		}

              }else if (evtNumHorDetSumGate>0){
		tigE_horDetSumGate->Fill(add_hit->GetEnergy());
		for (int horHit = 0; horHit < evtNumHorDetSumGate;horHit++){ //loop over all hits within the horizontal gate and plot the TIP channel number
		  tipChannel_HorizontalCut->Fill(tip_hit->GetTipChannel());
		}
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

          //TIGRESS-TIGRESS addback
          for (int tigHitIndAB2 = tigHitIndAB+1; tigHitIndAB2 < tigress->GetAddbackMultiplicity(); tigHitIndAB2++)
          {
            add_hit2 = tigress->GetAddbackHit(tigHitIndAB2);
            suppAdd = add_hit2->BGOFired();
            if (!suppAdd && add_hit2->GetEnergy() > 15)
            {
              addT_addT->Fill(add_hit->GetTime() - add_hit2->GetTime());
              addE_addE->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
              addE_addE->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
              if (gate1D((add_hit->GetTime() - add_hit2->GetTime()), tigtigT[0], tigtigT[1]))
              {
                addE_addE_tg->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                addE_addE_tg->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
              }
            }
          }
          
        }
      }
    }
    if(tip && tigress){
      numTipTigEvt++;
      
      if((tip->GetMultiplicity()>=0)||(tigress->GetAddbackMultiplicity()>=0)){
        Int_t tigSuppMult=0;
        for (int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
          add_hit = tigress->GetAddbackHit(tigHitIndAB);
          suppAdd = add_hit->BGOFired();
          if (!suppAdd && add_hit->GetEnergy() > 15){
            tigSuppMult++;
          }
        }
        tiptig_mult->Fill(tip->GetMultiplicity(),tigSuppMult);
      }
    }
    

    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree

  for(int i=0; i<NTIPRING; i++){
    tipPileupVsRing->SetBinContent(i+1,(float)(numTipRingPileupHits[i])/(float)(numTipRingHits[i]));
  }

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Number of TIGRESS addback hits: " << numTigressABHits << endl;
  cout << "Number of TIP hits: " << numTipHits << endl;
  cout << "Number of TIP pileup hits: " << numTipPileupHits << " (" << (float)(numTipPileupHits)/float(numTipPileupHits + numTipHits) << " %)" << endl;
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

  TDirectory *tipPIDGateddir = myfile->mkdir("TIP PID Gated");
  tipPIDGateddir->cd();
  tipPIDGatedList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  TDirectory *tigbgodir = myfile->mkdir("TIGRESS-BGO");
  tigbgodir->cd();
  tigbgoList->Write();
  myfile->cd();

  TDirectory *tiptigdir = myfile->mkdir("TIP-TIGRESS");
  tiptigdir->cd();
  tiptigList->Write();
  myfile->cd();

  myfile->Write();
  myfile->Close();
}
int main(int argc, char **argv)
{

  SortCode *mysort = new SortCode();

  char const *afile;
  char const *outfile;
  char const *calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0)
  {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1)
  {
    cout << "Arguments: SortData analysis_tree calibration_file output_file" << endl;
    cout << "Default values will be used if arguments (other than analysis tree) are omitted." << endl;
    return 0;
  }
  else if (argc == 2)
  {
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
  }
  else
  {
    printf("Too many arguments\nArguments: SortData analysis_tree calibration_file output_file\n");
    return 0;
  }

  theApp=new TApplication("App", &argc, argv);

  mysort->SortData(afile, calfile, outfile);

  return 0;
}
