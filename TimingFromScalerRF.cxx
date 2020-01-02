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
#include "TFragment.h"
#include "TEpicsFrag.h"
#include "TParserLibrary.h"
#include "TEnv.h"

#include "TEmma.h"
#include "TEmmaHit.h"
#include "TTigress.h"
#include "TTigressHit.h"
#include "TS3.h"
#include "TS3Hit.h"

#include <cstddef>

TH1D *periodHist, *phaseHist, *parHist;

std::vector<double>    phaseShiftBuffer, phaseBuffer, freqBuffer, parBuffer[4];
std::vector<ULong64_t> timestampBuffer;

ULong64_t mints, maxts;

char str[256];
char ftreename[256], atreename[256];

// builds a table of RF time values
// returns the number of RF scaler fragments in the table
int BuildRFTimeTable(TTree* stree)
{

   TEpicsFrag* currentFrag = nullptr;
   TBranch*    branch      = stree->GetBranch("TEpicsFrag");
   branch->SetAddress(&currentFrag);

   int entries = stree->GetEntries();

   int  FragsIn   = 0;
   int  RFFragsIn = 0;
   bool fe, p0e, p1e, p2e, de, te; // flags to check whether specific RF data exists in a fragment

   maxts = 0;
   mints = 1E16;

   // reserve enough space for long runs
   freqBuffer.reserve(entries);
   parBuffer[0].reserve(entries);
   parBuffer[1].reserve(entries);
   parBuffer[2].reserve(entries);
   parBuffer[3].reserve(entries);
   timestampBuffer.reserve(entries);

   for(int i = 0; i < entries; i++) {

      stree->GetEntry(i);
      FragsIn++;

      fe  = false;
      p0e = false;
      p1e = false;
      p2e = false;
      de  = false;
      double val;

      timestampBuffer.push_back(currentFrag->GetTimeStamp());

      if(currentFrag->GetTimeStamp() < mints) mints = currentFrag->GetTimeStamp();
      if(currentFrag->GetTimeStamp() > maxts) maxts = currentFrag->GetTimeStamp();

      for(int j = 0; j < currentFrag->fName.size(); j++) {
         val = currentFrag->GetData(j);
         // printf("%s: %f\n",currentFrag->fName[j].c_str(),val);
         if(j < currentFrag->fData.size()) {
            if(strcmp(currentFrag->fName[j].c_str(), "RF Frequency") == 0) {
               freqBuffer.push_back(val);
               fe = true;
            } else if(strcmp(currentFrag->fName[j].c_str(), "RF par 0") == 0) {
               parBuffer[0].push_back(val);
               p0e = true;
            } else if(strcmp(currentFrag->fName[j].c_str(), "RF par 1") == 0) {
               parBuffer[1].push_back(val);
               p1e = true;
            } else if(strcmp(currentFrag->fName[j].c_str(), "RF par 2") == 0) {
               parBuffer[2].push_back(val);
               p2e = true;
            } else if(strcmp(currentFrag->fName[j].c_str(), "RF Determinant") == 0) {
               parBuffer[3].push_back(val);
               de = true;
            }
         }
      }

      if(fe && p0e && p1e && p2e && de) {
         RFFragsIn++;
      } else {
         // make all vectors the same size and then remove the last element of each (corresponding to the bad fragment)
         if(!fe) {
            freqBuffer.push_back(0.0);
         }
         if(!p0e) {
            parBuffer[0].push_back(0.0);
         }
         if(!p1e) {
            parBuffer[1].push_back(0.0);
         }
         if(!p2e) {
            parBuffer[2].push_back(0.0);
         }
         if(!de) {
            parBuffer[3].push_back(0.0);
         }
         freqBuffer.pop_back();
         parBuffer[0].pop_back();
         parBuffer[1].pop_back();
         parBuffer[2].pop_back();
         parBuffer[3].pop_back();
         timestampBuffer.pop_back();
      }
   }

   printf("%i RF fragment(s) read in.\n", FragsIn);

   double t0, T, A, s, c;
   double par[3];

   double prevTs = 0.;
   double prevBadTs = 0.;
   double prevPhase = 0.;
   int previ = 0;

   for(int i = 0; i < RFFragsIn; i++) {

      if(parBuffer[3][i] != 0.) {
         par[0] = parBuffer[0][i] / parBuffer[3][i];
         par[1] = parBuffer[1][i] / parBuffer[3][i];
         par[2] = parBuffer[2][i] / parBuffer[3][i];
         T      = 1.34217728E9 / freqBuffer[i]; // period in ns
         // printf("parBuffer[0][i]: %f, parBuffer[1][i]: %f, parBuffer[2][i]: %f, parBuffer[3][i]:
         // %f\n",parBuffer[0][i],parBuffer[1][i],parBuffer[2][i],parBuffer[3][i]); printf("par0: %f, par1: %f, par2:
         // %f, T: %f\n",par[0],par[1],par[2],T); getc(stdin);
         A = sqrt(par[1] * par[1] + par[0] * par[0]);
         s = -par[0] / A;
         c = par[1] / A;
         if(s >= 0) {
            t0 = acos(c) * T / (2 * TMath::Pi());
         } else {
            t0 = (1 - acos(c) / (2 * TMath::Pi())) * T;
         }
         phaseShiftBuffer.push_back(t0); // in ns
         phaseBuffer.push_back((t0 / T) * (2 * TMath::Pi()));
         phaseHist->Fill(t0);
         periodHist->Fill(T); // period in ns
         parHist->Fill(A);
         //printf("i=%i, ts=%llu, phase=%f rad, phase shift=%f ns, T=%f ns\n",i, timestampBuffer[i], (t0 / T) * (2 * TMath::Pi()), t0, T);
         // getc(stdin);
         /*double tsms = timestampBuffer[i]/1000000.;
         double phaseDiff = phaseBuffer[i]-prevPhase;
         if(phaseDiff < 0)
            phaseDiff += (2 * TMath::Pi());
         hist->Fill(tsms,tsms-prevTs);
         hist2->Fill(tsms,phaseDiff);
         //printf("vals: %f %f\n",tsms,tsms-prevTs);
         printf("vals: %f %f\n",tsms,phaseDiff);
         if(phaseDiff< 4.0){
            hist3->Fill(tsms-prevBadTs);
            prevBadTs=tsms;
         }  
         //getc(stdin);
         prevTs=tsms;
         
         prevPhase=phaseBuffer[i];*/

      }else{
         printf("Skipping failed RF fit.\n"); 
      }

      
      
   }

   // getc(stdin);
   return RFFragsIn;
}

void PhaseInterpTest(int numRFFrags)
{

   TH1D* interpHist =
      new TH1D("interpolated - actual phase shift", "interpolated - actual phase shift", 20000, -100, 100);
   interpHist->Reset();

   double phase, interpPhaseShift, interpPhase;
   double T, numTinterp;

   for(int i = 1; i < numRFFrags; i++) {
      phase = phaseBuffer[i];
      T     = 1.34217728E9 / freqBuffer[i]; // period in ns
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

double GetPhaseNsForTimestamp(ULong64_t timestamp, int numRFFrags)
{
   // printf("ts: %llu\n",timestamp);
   double T, numTinterp, interpPhaseShift;
   double interpPhase = -100;
   int    i           = 0;
   for(i = 0; i < numRFFrags; i++) {
      if((i > 0) && (timestampBuffer[i] > timestamp) &&
         (timestampBuffer[i - 1] <= timestamp)) { // ignoring rollover for now
         T = 1.34217728E9 / freqBuffer[i - 1];    // period in ns
         if(timestampBuffer[i] > timestampBuffer[i - 1])
            numTinterp = (timestamp - timestampBuffer[i - 1]) / T;
         else // handle rollover of timestamps
            numTinterp =
               ((timestamp + 43980465111040) - timestampBuffer[i - 1]) / T; // timestamp is a 42-bit unsigned int
         interpPhaseShift = (numTinterp - (int)numTinterp) * T;
         if(phaseShiftBuffer[i - 1] + interpPhaseShift < T)
            interpPhase = phaseShiftBuffer[i - 1] + interpPhaseShift;
         else
            interpPhase = interpPhaseShift - (T - phaseShiftBuffer[i - 1]);
         break;
      }
   }
   // if((interpPhase > 66.16) && (interpPhase < 66.175)) printf("ok ts: %llu, phase: %f\n", timestamp, interpPhase);
   // printf("numTinterp: %f, i: %i, interpPhaseShift: %f, interpPhase: %f\n",numTinterp, i, interpPhaseShift,
   // interpPhase);
   return interpPhase;
}

// below: if true, interpolate from below, if false, interpolate from above
double GetPhaseRadForTimestamp(ULong64_t timestamp, int numRFFrags, bool below)
{
   // printf("ts: %llu\n",timestamp);
   double T;
   double numTinterp, interpPhaseShift;
   double interpPhase = -100;
   int    i           = 0;
   for(i = 0; i < numRFFrags; i++) {
      if((i > 0) && (timestampBuffer[i] > timestamp) &&
         (timestampBuffer[i - 1] <= timestamp)) { // ignoring rollover for now
         if(below) { //interpolate from previous RF event
            T                = 1.34217728E9 / (freqBuffer[i - 1]); // period in ns
            numTinterp       = (timestamp - timestampBuffer[i - 1]) / T; //number of periods to interpolate
            interpPhaseShift = (numTinterp - (int)numTinterp) * (2 * TMath::Pi()); //phase shift needed
            // set phase, wrapping around so that it ranges from 0 to 2pi
            if((phaseBuffer[i - 1] + interpPhaseShift) < (2 * TMath::Pi()))
               interpPhase = phaseBuffer[i - 1] + interpPhaseShift;
            else
               interpPhase = phaseBuffer[i - 1] + interpPhaseShift - (2 * TMath::Pi());
            break;
         } else { //interpolate from RF event after
            T                = 1.34217728E9 / freqBuffer[i]; // period in ns
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
   // printf("numTinterp: %f, i: %i, interpPhaseShift: %f, interpPhase: %f\n",numTinterp, i, interpPhaseShift,
   // interpPhase);
   return interpPhase;
}

void MapPhaseTest(TTree* ftree, int numRFFrags, TH1D* interpHist, int chan1, int chan2, bool append)
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
         double phase1 = GetPhaseRadForTimestamp(ts1, numRFFrags, true);
         // interpHist->Fill(fmod(correctedTS,84.841652));
         interpHist->Fill(phase1);
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

void MapPhaseTest2D(TTree* ftree, int numRFFrags, TH2D* interpHist, int chan1, int chan2, bool append)
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
         double phase1 = GetPhaseRadForTimestamp(ts1, numRFFrags, true);

         for(int j = i; j < (i + 3000); j++) {
            ftree->GetEntry(j);
            if(currentFrag->GetChannelNumber() == chan2) {
               ts2           = currentFrag->GetTime();
               double phase2 = GetPhaseRadForTimestamp(ts2, numRFFrags, true);
               interpHist->Fill(phase1, phase2);
            }
         }
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

void MapPhaseProgressionTest(TTree* atree, int numRFFrags, TH2D* interpHist, int hitType, bool append)
{

   if(!append) interpHist->Reset();

   ULong64_t ts;
   int entries = atree->GetEntries();
   //int entries = 200000;
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
         for (int j = 0; j < emma->GetSSBMultiplicity(); j++) { // Get SSB hits
            ssb_hit = emma->GetSSBHit(j);
            ts = ssb_hit->GetTime();
            double phase = GetPhaseRadForTimestamp(ts, numRFFrags, true); //get RF phase
            interpHist->Fill(phase, ssb_hit->GetEnergy());
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
         for (int j = 0; j < tig->GetMultiplicity(); j++) { // Get TIGRESS hits
            tig_hit = tig->GetTigressHit(j);
            ts = tig_hit->GetTime();
            double phase = GetPhaseRadForTimestamp(ts, numRFFrags, true); //get RF phase
            interpHist->Fill(phase, tig_hit->GetEnergy());
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
         for (int j = 0; j < s3->GetPixelMultiplicity(); j++) { // Get S3 hits
            s3_hit = s3->GetPixelHit(j);
            ts = s3_hit->GetTime();
            double phase = GetPhaseRadForTimestamp(ts, numRFFrags, true); //get RF phase
            interpHist->Fill(phase, s3_hit->GetEnergy());
         }
      }
   }else{
      printf("Invalid hit type!\n");
      exit(-1);
   }

   

   

   // getc(stdin);
}

void PhaseConsistencyTest(TTree* ftree, int numRFFrags, TH1D* interpHist, int chan1, int chan2, bool append)
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
            GetPhaseRadForTimestamp(ts1, numRFFrags, true) - GetPhaseRadForTimestamp(ts1, numRFFrags, false);
         // interpHist->Fill(fmod(correctedTS,84.841652));
         interpHist->Fill(phasediff);
      }

      // currentFrag->Print();
   }

   // getc(stdin);
}

/* void MapPhase2DTest(TTree* atree, int numRFFrags, TH2D *hist, bool append){

   if(!append)
      hist->Reset();

   ULong64_t ts;
   //double modts1,modts2;
   double tigE;

   //TEmma* currentEmmaFrag = nullptr;
   TTigress* currentTigFrag = nullptr;
   //TBranch* branchEmma = atree->GetBranch("TEmma");
   //branchEmma->SetAddress(&currentEmmaFrag);
   TBranch* branchTig = atree->GetBranch("TTigress");
   branchTig->SetAddress(&currentTigFrag);

   int entries = atree->GetEntries();

   for(int i=0;i<entries;i++){

      if(i%100==0)
         printf("Entry %i / %i\r",i,entries);

      atree->GetEntry(i);

      for(int j=0;j<currentTigFrag->GetMultiplicity();j++){
         tigE = currentTigFrag->GetTigressHit(j)->GetEnergy()/1000.0;
         ts = static_cast<ULong64_t>(currentTigFrag->GetTigressHit(j)->GetTimeStampNs());
         hist->Fill(tigE,GetPhaseNsForTimestamp(ts,numRFFrags));
      }

   }

} */

/* void WriteCorrectedTimes(TTree* inptree, TTree* outtree)
{

   TBranch* branch = inptree->GetBranch("TTigress");
   if(branch == nullptr) {
      printf("Tigress branch not in analysis tree.\n");
   } else {
      printf("Tigress branch in analysis tree.\n");
      TTigress* currentDet = nullptr;
      branch->SetAddress(&currentDet);

      int entries = inptree->GetEntries();

      int FragsIn = 0;

      for(int i=0;i<entries;i++){

         inptree->GetEntry(i);
         FragsIn++;
         printf("Fragment %i\n",FragsIn);
         for(int j=0;j<currentDet->GetHitVector().size();j++){
            printf("ts=%li\n",currentDet->GetHitVector()[j]->GetTimeStampNs());
         }
         getc(stdin);

      }

   }

   return;
} */

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
      printf("\nThis code plots the RF phase for a given hit type in the specified run.\nRF data should be present in the fragment epics tree(s). Valid hit types: s3, ssb, tigress\n");
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
   TH2D* interpHist2 = new TH2D("Timing", "Timing", 70, (-1 * TMath::Pi()), (1 * TMath::Pi()), 1000, 0,70000);
   interpHist2->Reset();
   /*TH2D* tsdiffvst = new TH2D("Timestamp difference vs timestamp", "Timestamp difference vs timestamp", 10000, 0, 1000000, 200, 0, 20);
   tsdiffvst->Reset();
   TH2D* phasediffvst = new TH2D("RF Phase difference vs timestamp", "RF Phase difference vs timestamp", 10000, 0, 1000000, 200, (0 * TMath::Pi()), (2 * TMath::Pi()));
   phasediffvst->Reset();*/

   int runNum = atoi(argv[1]);

   int numTrees = 0;
   bool searching = true;
   while(searching) {
      sprintf(ftreename, "fragment%i_%03d.root",runNum,numTrees);
      sprintf(atreename, "analysis%i_%03d.root",runNum,numTrees);
      ffile = new TFile(ftreename);
      afile = new TFile(atreename);
      if((ffile == nullptr) || (!ffile->IsOpen())) {
         printf("%i tree(s) found for run %i.\n", numTrees,runNum);
         searching = false;
      }else if((afile == nullptr) || (!afile->IsOpen())) {
         printf("%i tree(s) found for run %i.\n", numTrees,runNum);
         searching = false;
      }else{
         printf("Files %s and %s found.\n",ftreename,atreename);
         numTrees++;
      }
      ffile->Close();
   }

   int hitType = -1;
   if(strcmp(argv[3],"ssb")==0){
      hitType = 1;
      printf("Will plot timing for SSB (from EMMA) hits.\n");
   }else if(strcmp(argv[3],"tigress")==0){
      hitType = 2;
      printf("Will plot timing for TIGRESS hits.\n");
   }else if(strcmp(argv[3],"s3")==0){
      hitType = 3;
      printf("Will plot timing for S3 hits.\n");
   }else{
      printf("ERROR: unknown hit type specified. Options are:\ns3 - S3 pixel hits\nssb - SSB (from EMMA)\ntigress - TIGRESS cores\n");
      exit(-1);
   }

   // scan the list files for ROOT files
   for(int i = 0; i < numTrees; i++) {

      sprintf(ftreename, "fragment%i_%03d.root",runNum,i);
      sprintf(atreename, "analysis%i_%03d.root",runNum,i);

      ffile = new TFile(ftreename);
      if((ffile == nullptr) || (!ffile->IsOpen())) {
         printf("Failed to open file '%s'!\n", ftreename);
         exit(-1);
      }

      TTree* epicstree = dynamic_cast<TTree*>(ffile->Get("EpicsTree"));
      if(epicstree == nullptr) {
         printf("Failed to find epics tree in file: %s\n", ftreename);
         exit(-1);
      } else {
         std::cout << epicstree->GetEntries() << " epics tree entries" << std::endl;
      }

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
      phaseShiftBuffer = {};
      freqBuffer       = {};
      for(int i = 0; i < 4; i++) parBuffer[i] = {};

      periodHist = new TH1D("period", "period", 20000, 84.83, 84.85);
      periodHist->Reset();
      phaseHist = new TH1D("phase", "phase", 20000, 0, 200);
      phaseHist->Reset();
      parHist = new TH1D("parameter 0 1 sumsquare", "parameter 0 1 sumsquare", 20000, -20, 20);
      parHist->Reset();

      int numRFFrags = BuildRFTimeTable(epicstree);
      //PhaseInterpTest(numRFFrags);

      // PhaseConsistencyTest(inptree, numRFFrags, interpHist, chan1, chan2, true);
      // MapPhaseTest(inptree, numRFFrags, interpHist, chan1, chan2, true);
      // MapPhaseTest2D(inptree, numRFFrags, interpHist2, chan1, chan2, true);
      MapPhaseProgressionTest(inptree, numRFFrags, interpHist2, hitType, true);
      // MapPhase2DTest(inptree,numRFFrags,hist,true);

      ffile->Close();
      afile->Close();
   }

   // TCanvas* c1 = new TCanvas("c1","c1",800,600);
   // periodHist->Draw();
   // phaseHist->Draw();
   // parHist->Draw();
   theApp      = new TApplication("App", &argc, argv);
   TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
   //interpHist->GetXaxis()->SetTitle("ts diff(ms)");
   //interpHist->Draw();
   sprintf(str,"%s RF phase vs. E",argv[3]);
   interpHist2->SetTitle(str);
   interpHist2->GetXaxis()->SetTitle("phase (rad)");
   interpHist2->GetYaxis()->SetTitle("Energy (keV)");
   interpHist2->Draw("colz");
   /*tsdiffvst->GetXaxis()->SetTitle("timestamp (ms)");
   tsdiffvst->GetYaxis()->SetTitle("timestamp diff (ms)");*/
   //tsdiffvst->Draw("colz");
   /*phasediffvst->GetXaxis()->SetTitle("timestamp (ms)");
   phasediffvst->GetYaxis()->SetTitle("phase diff (rad)");
   phasediffvst->Draw("colz");*/
   theApp->Run(kTRUE);

   // WriteCorrectedTimes(inptree,outtree);

   // auto* outfile = new TFile(argv[3], "recreate");
   // list->Write();

   return 0;
}

#endif
