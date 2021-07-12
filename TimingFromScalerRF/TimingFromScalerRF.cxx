// g++ TimingFromScalerRF.cxx -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -o TimingFromScalerRF

// Some test functions for dealing with RF data
//
//
#include <utility>
#include <vector>
#include <cstdio>
#include <iostream>
#include <iomanip>

#include "TTree.h"
#include "TTreeIndex.h"
#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"

#include "TChannel.h"
#include "TDetectorHit.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TRF.h"
#include "TEmma.h"
#include "TEmmaHit.h"
#include "TTigress.h"
#include "TTigressHit.h"
#include "TS3.h"
#include "TS3Hit.h"

#include <cstddef>

TH1D *periodHist, *phaseHist, *parHist;

std::vector<double>    phaseBuffer, periodBuffer;
std::vector<double>   timestampBuffer;

char str[256];
char atreename[256];

// builds a table of RF time values
// returns the number of RF scaler fragments in the table
int BuildRFTimeTable(TTree* atree)
{

   TRF* rf = nullptr;
   TBranch*    branch      = atree->GetBranch("TRF");
   branch->SetAddress(&rf);

   int entries = atree->GetEntries();

   int  EvtsIn   = 0;

   // reserve enough space for long runs
   timestampBuffer.reserve(entries);
   periodBuffer.reserve(entries);
   phaseBuffer.reserve(entries);

   for(int i = 0; i < entries; i++) {

      atree->GetEntry(i);
      EvtsIn++;

      //std::cout << "ts: " << rf->TimeStamp() << ", per: " << rf->Period() << ", phase: " << rf->Phase() << std::endl;
      if((rf->TimeStamp() > 0)&&(rf->Phase() >= 0.0)){
         timestampBuffer.push_back(rf->TimeStamp() * 10.0);
         periodBuffer.push_back(rf->Period());
         phaseBuffer.push_back(rf->Phase());
      }
      
   }

   printf("%i RF event(s) read in.\n", EvtsIn);

   // getc(stdin);
   return EvtsIn;
}

void PhaseInterpTest(int numRFEvt)
{

   TH1D* interpHist =
      new TH1D("interpolated - actual phase shift", "interpolated - actual phase shift", 20000, -100, 100);
   interpHist->Reset();

   double phase, interpPhaseShift, interpPhase;
   double T, numTinterp;

   for(int i = 1; i < numRFEvt; i++) {
      phase = phaseBuffer[i];
      T     = periodBuffer[i]; // period in ns
      if(timestampBuffer[i] > timestampBuffer[i - 1])
         numTinterp = (timestampBuffer[i] - timestampBuffer[i - 1]) / T;
      else // handle rollover of timestamps
         numTinterp =
            ((timestampBuffer[i] + 43980465111040) - timestampBuffer[i - 1]) / T; // timestamp is a 42-bit unsigned int
      interpPhaseShift = (numTinterp - (int)numTinterp) * (2 * TMath::Pi());
      if(phaseBuffer[i - 1] + interpPhaseShift < (2 * TMath::Pi()))
         interpPhase = phaseBuffer[i - 1] + interpPhaseShift;
      else
         interpPhase = interpPhaseShift - ((2 * TMath::Pi()) - phaseBuffer[i - 1]);
       printf("phase: %f, numTinterp: %f, interpPhase: %f\n",phase,numTinterp,interpPhase);
       //getc(stdin);
      interpHist->Fill(interpPhase - phase);
   }

   TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
   interpHist->GetXaxis()->SetTitle("interpolated - actual phase shift (ns)");
   interpHist->Draw();
}

// below: if true, interpolate from below, if false, interpolate from above
double GetPhaseRadForTimestamp(ULong64_t timestamp, int numRFEvt, bool below)
{
   // printf("ts: %llu\n",timestamp);
   double T;
   double numTinterp, interpPhaseShift;
   double interpPhase = -100;
   int    i           = 0;
   for(i = 0; i < numRFEvt; i++) {
      if((i > 0) && (timestampBuffer[i] > timestamp) &&
         (timestampBuffer[i - 1] <= timestamp)) { // ignoring rollover for now
         if(below) { //interpolate from previous RF event
            T                = periodBuffer[i - 1]; // period in ns
            numTinterp       = (timestamp - timestampBuffer[i - 1]) / T; //number of periods to interpolate
            interpPhaseShift = (numTinterp - (int)numTinterp) * (2 * TMath::Pi()); //phase shift needed
            // set phase, wrapping around so that it ranges from 0 to 2pi
            if((phaseBuffer[i - 1] + interpPhaseShift) < (2 * TMath::Pi()))
               interpPhase = phaseBuffer[i - 1] + interpPhaseShift;
            else
               interpPhase = phaseBuffer[i - 1] + interpPhaseShift - (2 * TMath::Pi());
            break;
         } else { //interpolate from RF event after
            T                = periodBuffer[i]; // period in ns
            numTinterp       = (timestampBuffer[i] - timestamp) / T; //number of periods to interpolate
            interpPhaseShift = (numTinterp - (int)numTinterp) * (2 * TMath::Pi()); //phase shift needed
            // set phase, wrapping around so that it ranges from 0 to 2pi
            if((phaseBuffer[i] - interpPhaseShift) < (2 * TMath::Pi()))
               interpPhase = phaseBuffer[i] - interpPhaseShift;
            else
               interpPhase = phaseBuffer[i] - interpPhaseShift - (2 * TMath::Pi());
            break;
         }
      }
   }
   // re-center phase to range from -pi to pi rather than 0 to 2pi
   if(interpPhase>TMath::Pi()){
       interpPhase-=2*TMath::Pi();
   }
   //printf("numTinterp: %f, i: %i, interpPhaseShift: %f, interpPhase: %f\n",numTinterp, i, interpPhaseShift, interpPhase);
   return interpPhase;
}

void MapPhaseTest(TTree* ftree, int numRFEvt, TH1D* interpHist, int chan1, int chan2, bool append)
{

   if(!append) interpHist->Reset();

   ULong64_t ts1, ts2;

   TFragment* currentFrag = nullptr;
   // TDetectorHit* currentHit = nullptr;
   TBranch* fragBranch = ftree->GetBranch("TFragment");
   // TBranch* hitBranch = ftree->GetBranch("TDetectorHit");
   fragBranch->SetAddress(&currentFrag);
   // hitBranch->SetAddress(&currentHit);

   int entries = ftree->GetEntries();

   int FragsIn = 0;

   for(int i = 30550000; i < 30850000; i++) {

      if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
      // printf("Entry %i / %i\nts: %llu\n",i,entries,currentFrag->GetTimeStampNs());

      ftree->GetEntry(i);
      FragsIn++;

      if(currentFrag->GetChannelNumber() == chan1) {
         // printf("channel address: %u\n",currentFrag->GetChannelNumber());
         ts1 = currentFrag->GetTime();
         // printf("ts: %llu\n",ts);
         double phase1 = GetPhaseRadForTimestamp(ts1, numRFEvt, true);
         // interpHist->Fill(fmod(correctedTS,84.841652));
         interpHist->Fill(phase1);
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

void MapPhaseTest2D(TTree* ftree, int numRFEvt, TH2D* interpHist, int chan1, int chan2, bool append)
{

   if(!append) interpHist->Reset();

   ULong64_t ts1, ts2;

   TFragment* currentFrag = nullptr;
   // TDetectorHit* currentHit = nullptr;
   TBranch* fragBranch = ftree->GetBranch("TFragment");
   // TBranch* hitBranch = ftree->GetBranch("TDetectorHit");
   fragBranch->SetAddress(&currentFrag);
   // hitBranch->SetAddress(&currentHit);

   int entries = ftree->GetEntries();

   int FragsIn = 0;

   for(int i = 50000; i < 1600000; i++) {

      if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
      // printf("Entry %i / %i\nts: %llu\n",i,entries,currentFrag->GetTimeStampNs());

      ftree->GetEntry(i);
      FragsIn++;

      if(currentFrag->GetChannelNumber() == chan1) {
         // printf("channel address: %u\n",currentFrag->GetChannelNumber());
         ts1 = currentFrag->GetTime();
         // printf("ts: %llu\n",ts);
         double phase1 = GetPhaseRadForTimestamp(ts1, numRFEvt, true);

         for(int j = i; j < (i + 3000); j++) {
            ftree->GetEntry(j);
            if(currentFrag->GetChannelNumber() == chan2) {
               ts2           = currentFrag->GetTime();
               double phase2 = GetPhaseRadForTimestamp(ts2, numRFEvt, true);
               interpHist->Fill(phase1, phase2);
            }
         }
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

void MapPhaseEnergyTest(TTree* atree, int numRFEvt, TH2D* interpHist, int hitType, bool append)
{

   if(!append) interpHist->Reset();

   ULong64_t ts;
   //int entries = atree->GetEntries();
   int entries = 600000;
   if(entries > atree->GetEntries())
     entries = atree->GetEntries();
   int entriesIn = 0;

   if(hitType==1){
      TEmma* emma = nullptr;
      TEmmaHit *ssb_hit;
      TBranch* anBranch = atree->GetBranch("TEmma");
      anBranch->SetAddress(&emma);
     
      for(int i = 0; i < entries; i++) {
         if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
         atree->GetEntry(i);
         entriesIn++;
         if(emma){
            for (int j = 0; j < emma->GetSSBMultiplicity(); j++) { // Get SSB hits
               ssb_hit = emma->GetSSBHit(j);
               ts = ssb_hit->GetTime();
               double phase = GetPhaseRadForTimestamp(ts, numRFEvt, true); //get RF phase
               interpHist->Fill(phase, ssb_hit->GetEnergy());
            }
         }
      }
   }else if (hitType==2){
      TTigress* tig = nullptr;
      TTigressHit *tig_hit;
      TBranch* anBranch = atree->GetBranch("TTigress");
      anBranch->SetAddress(&tig);

      for(int i = 0; i < entries; i++) {
         if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
         atree->GetEntry(i);
         entriesIn++;
         if(tig){
            for (int j = 0; j < tig->GetMultiplicity(); j++) { // Get TIGRESS hits
               tig_hit = tig->GetTigressHit(j);
               ts = tig_hit->GetTime();
               double phase = GetPhaseRadForTimestamp(ts, numRFEvt, true); //get RF phase
               interpHist->Fill(phase, tig_hit->GetEnergy());
            }
         }
      }
   }else if (hitType==3){
      TS3* s3 = nullptr;
      TS3Hit *s3_hit;
      TBranch* anBranch = atree->GetBranch("TS3");
      anBranch->SetAddress(&s3);

      for(int i = 0; i < entries; i++) {
         if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
         atree->GetEntry(i);
         entriesIn++;
         if(s3){
            for (int j = 0; j < s3->GetPixelMultiplicity(); j++) { // Get S3 hits
               s3_hit = s3->GetPixelHit(j);
               ts = s3_hit->GetTime();
               double phase = GetPhaseRadForTimestamp(ts, numRFEvt, true); //get RF phase
               interpHist->Fill(phase, s3_hit->GetEnergy());
            }
         }
      }
   }else if(hitType==4){
      TEmma* emma = nullptr;
      TEmmaHit *ssb_hit;
      TBranch* anBranch = atree->GetBranch("TEmma");
      anBranch->SetAddress(&emma);
      TTigress* tig = nullptr;
      TTigressHit *tig_hit;
      TBranch* anBranch2 = atree->GetBranch("TTigress");
      anBranch2->SetAddress(&tig);
     
      for(int i = 0; i < entries; i++) {
         if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
         atree->GetEntry(i);
         entriesIn++;
         if(emma&&tig){
            for (int j = 0; j < emma->GetSSBMultiplicity(); j++) { // Get SSB hits
               ssb_hit = emma->GetSSBHit(j);
               ts = ssb_hit->GetTime();
               double phase = GetPhaseRadForTimestamp(ts, numRFEvt, true); //get RF phase
               interpHist->Fill(phase, ssb_hit->GetEnergy());
            }
         }
      }
   }else{
      printf("Invalid hit type!\n");
      exit(-1);
   }

   

   

   // getc(stdin);
}

void PhaseConsistencyTest(TTree* ftree, int numRFEvt, TH1D* interpHist, int chan1, int chan2, bool append)
{

   if(!append) interpHist->Reset();

   ULong64_t ts1, ts2;

   TFragment* currentFrag = nullptr;
   // TDetectorHit* currentHit = nullptr;
   TBranch* fragBranch = ftree->GetBranch("TFragment");
   // TBranch* hitBranch = ftree->GetBranch("TDetectorHit");
   fragBranch->SetAddress(&currentFrag);
   // hitBranch->SetAddress(&currentHit);

   int entries = ftree->GetEntries();

   int FragsIn = 0;

   for(int i = 0; i < entries; i++) {

      if(i % 100 == 0) printf("Entry %i / %i\r", i, entries);
      // printf("Entry %i / %i\nts: %llu\n",i,entries,currentFrag->GetTimeStampNs());

      ftree->GetEntry(i);
      FragsIn++;

      if(currentFrag->GetChannelNumber() == chan1) {
         // printf("channel address: %u\n",currentFrag->GetChannelNumber());
         ts1 = currentFrag->GetTime();
         // printf("ts: %llu\n",ts);
         double phasediff =
            GetPhaseRadForTimestamp(ts1, numRFEvt, true) - GetPhaseRadForTimestamp(ts1, numRFEvt, false);
         // interpHist->Fill(fmod(correctedTS,84.841652));
         interpHist->Fill(phasediff);
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

#ifndef __CINT__

int main(int argc, char** argv)
{
   TApplication* theApp;

   FILE*  list;
   TFile *ffile, *afile;
   TFile* inpfile;
   char const * calfile;

   if(argc != 4) {
      printf("%s run_number calfile hit_type\n", argv[0]);
      printf("\nThis code plots the RF phase for a given hit type in the specified run.\nRF data should be present in the analysis tree(s). Valid hit types: s3, ssb, tigress, ssb_tigress (SSB in coinc with TIGRESS)\n");
      return 0;
   }

  //load configuration, needed to access channel numbers, and to load the proper time stamp units
  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0) {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

   calfile = argv[2]; 
   printf("Reading calibration file: %s\n", calfile);
   TChannel::ReadCalFile(calfile);


   TH1D* interpHist = new TH1D("Timing", "Timing", 5000, 0, 2048);
   interpHist->Reset();
   TH2D* interpHist2;
   
   /*TH2D* tsdiffvst = new TH2D("Timestamp difference vs timestamp", "Timestamp difference vs timestamp", 10000, 0, 1000000, 200, 0, 20);
   tsdiffvst->Reset();
   TH2D* phasediffvst = new TH2D("RF Phase difference vs timestamp", "RF Phase difference vs timestamp", 10000, 0, 1000000, 200, (0 * TMath::Pi()), (2 * TMath::Pi()));
   phasediffvst->Reset();*/

   int runNum = atoi(argv[1]);

   int numTrees = 0;
   bool searching = true;
   while(searching){
      sprintf(atreename, "analysis%i_%03d.root",runNum,numTrees);
      afile = new TFile(atreename);
      if((afile == nullptr) || (!afile->IsOpen())) {
         printf("%i tree(s) found for run %i.\n", numTrees,runNum);
         searching = false;
      }else{
         printf("File %s found.\n",atreename);
         numTrees++;
      }
      afile->Close();
   }

   int hitType = -1;
   if(strcmp(argv[3],"ssb")==0){
      hitType = 1;
      interpHist2 = new TH2D("Timing", "Timing", 100, (-1 * TMath::Pi()), (1 * TMath::Pi()), 1000, 0,80000);
      printf("Will plot timing for SSB (from EMMA) hits.\n");
   }else if(strcmp(argv[3],"tigress")==0){
      hitType = 2;
      interpHist2 = new TH2D("Timing", "Timing", 100, (-1 * TMath::Pi()), (1 * TMath::Pi()), 1000, 0,8000);
      printf("Will plot timing for TIGRESS hits.\n");
   }else if(strcmp(argv[3],"s3")==0){
      hitType = 3;
      interpHist2 = new TH2D("Timing", "Timing", 100, (-1 * TMath::Pi()), (1 * TMath::Pi()), 1000, 0,14000);
      printf("Will plot timing for S3 hits.\n");
   }else if(strcmp(argv[3],"ssb_tigress")==0){
      hitType = 4;
      interpHist2 = new TH2D("Timing", "Timing", 100, (-1 * TMath::Pi()), (1 * TMath::Pi()), 1000, 0,80000);
      printf("Will plot timing for SSB hits in coinc with TIGRESS hits.\n");
   }else{
      printf("ERROR: unknown hit type specified. Options are:\ns3 - S3 pixel hits\nssb - SSB (from EMMA)\ntigress - TIGRESS cores\n");
      exit(-1);
   }
   interpHist2->Reset();

   // scan the list files for ROOT files
   for(int i = 0; i < numTrees; i++) {

      sprintf(atreename, "analysis%i_%03d.root",runNum,i);

      afile = new TFile(atreename);
      if((afile == nullptr) || (!afile->IsOpen())) {
         printf("Failed to open file '%s'!\n", atreename);
         exit(-1);
      }

      TTree* inptree = dynamic_cast<TTree*>(afile->Get("AnalysisTree"));
      if(inptree == nullptr) {
         printf("Failed to find analysis tree in file: %s\n", atreename);
         exit(-1);
      } else {
         std::cout << inptree->GetEntries() << " analysis tree entries" << std::endl;
      }


      // initialize arrays
      timestampBuffer  = {};
      phaseBuffer      = {};
      periodBuffer     = {};

      periodHist = new TH1D("period", "period", 20000, 84.83, 84.85);
      periodHist->Reset();
      phaseHist = new TH1D("phase", "phase", 20000, 0, 200);
      phaseHist->Reset();
      parHist = new TH1D("parameter 0 1 sumsquare", "parameter 0 1 sumsquare", 20000, -20, 20);
      parHist->Reset();

      int numRFEvt = BuildRFTimeTable(inptree);
      //PhaseInterpTest(numRFEvt);

      // PhaseConsistencyTest(inptree, numRFEvt, interpHist, chan1, chan2, true);
      // MapPhaseTest(inptree, numRFEvt, interpHist, chan1, chan2, true);
      // MapPhaseTest2D(inptree, numRFEvt, interpHist2, chan1, chan2, true);
      MapPhaseEnergyTest(inptree, numRFEvt, interpHist2, hitType, true);
      // MapPhase2DTest(inptree,numRFEvt,hist,true);

      afile->Close();
   }

   // TCanvas* c1 = new TCanvas("c1","c1",800,600);
   // periodHist->Draw();
   // phaseHist->Draw();
   // parHist->Draw();
   theApp      = new TApplication("App", &argc, argv);
   TCanvas* c1 = new TCanvas("c1", "c1", 1000, 900);
   //interpHist->GetXaxis()->SetTitle("ts diff(ms)");
   //interpHist->Draw();
   sprintf(str,"%s RF phase vs. E",argv[3]);
   interpHist2->SetTitle(str);
   interpHist2->GetXaxis()->SetTitle("phase (rad)");
   interpHist2->GetYaxis()->SetTitle("S3 energy");
   interpHist2->Draw("colz");
   /*tsdiffvst->GetXaxis()->SetTitle("timestamp (ms)");
   tsdiffvst->GetYaxis()->SetTitle("timestamp diff (ms)");*/
   //tsdiffvst->Draw("colz");
   /*phasediffvst->GetXaxis()->SetTitle("timestamp (ms)");
   phasediffvst->GetYaxis()->SetTitle("phase diff (rad)");
   phasediffvst->Draw("colz");*/
   theApp->Run(kTRUE);

   // auto* outfile = new TFile(argv[3], "recreate");
   // list->Write();

   return 0;
}

#endif
