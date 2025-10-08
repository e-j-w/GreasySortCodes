//use the Makefile!

#define SortDiagnosticsS_cxx
#include "common.cxx"
#include "SortDiagnosticsSMOL.h"

using namespace std;


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
  uint64_t pileupCtrs[16];
  fread(&sentries,sizeof(uint64_t),1,inp);
  uint64_t smolVersion = (uint64_t)(sentries >> 48);
  if(smolVersion > 0){
    fread(&pileupCtrs,sizeof(pileupCtrs),1,inp);
    printf("\nNumber of hits of each pileup type:\n");
    uint64_t totalHits = 0;
    for(uint8_t i=0; i<16; i++){
      printf("Pileup type %2u: %Lu\n",i,pileupCtrs[i]);
      totalHits += pileupCtrs[i];
    }
    printf("Total hits:     %Lu\n",totalHits);
    long double frac = (long double)(pileupCtrs[1])/((long double)(totalHits));
    printf("Fraction of hits with type 1 (no pileup): %Lf\n",frac);
  }
  sentries &= 0xFFFFFFFFFFFF; // only first 48 bits specify number of events
  sorted_evt sortedEvt;
  
  uint8_t footerVal;

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

    /*if(jentry < 1008000000){
      if (jentry % 9713 == 0)
        cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
      continue;
    }*/

    /*if(jentry >= 5052150){
      printf("\n\n\n  event %Lu\n  metadata: %u\n",jentry,sortedEvt.header.metadata);
			printf("  event time: %f\n",sortedEvt.header.evtTimeNs);
			printf("  num hits: %u\n",sortedEvt.header.numNoABHits);
			for(int i = 0; i<sortedEvt.header.numNoABHits;i++){
				printf("    hit %u - time offset: %f, energy: %f, core: %u\n",i,sortedEvt.noABHit[i].timeOffsetNs,sortedEvt.noABHit[i].energy,sortedEvt.noABHit[i].core & 63U);
			}
    }*/

    hpgeMult->Fill(sortedEvt.header.numNoABHits);
    
    for (int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){
        hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy);
        hpgeE_ANum->Fill(sortedEvt.noABHit[noABHitInd].core & 63U, sortedEvt.noABHit[noABHitInd].energy);
      }
    }

    //check for events with gate conditions
    int gate1477 = 0;
    int gate1477HitInd = 0;
    int gate685 = 0;
    int gate685HitInd = 0;
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if((sortedEvt.noABHit[noABHitInd].energy >= 1470)&&(sortedEvt.noABHit[noABHitInd].energy < 1484)){
        gate1477 = 1;
        gate1477HitInd = noABHitInd;
        break;
      }
    }
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){
      if((sortedEvt.noABHit[noABHitInd].energy >= 680)&&(sortedEvt.noABHit[noABHitInd].energy < 690)){
        gate685 = 1;
        gate685HitInd = noABHitInd;
        break;
      }
    }

    //evaluate non-addback timing conditions
    uint8_t sameEvtCtr = 0;
    for(int noABHitInd = 0; noABHitInd < sortedEvt.header.numNoABHits; noABHitInd++){

      if(sortedEvt.noABHit[noABHitInd].energy > MIN_HPGE_EAB){

        //TIG-TIG timing, position, and energy
        for(int noABHitInd2 = noABHitInd+1; noABHitInd2 < sortedEvt.header.numNoABHits; noABHitInd2++){
          if(sortedEvt.noABHit[noABHitInd2].energy > MIN_HPGE_EAB){
            Double_t tDiff = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,noABHitInd);
            /*if(tDiff >= 0.0f && tDiff <= 1.0f){
              sameEvtCtr++;
              if(sameEvtCtr > 2){
                printf("tDiff: %f, event %li hits %i and %i\n",tDiff,jentry,noABHitInd,noABHitInd2);
                for(int noABHitInd3 = 0; noABHitInd3 < sortedEvt.header.numNoABHits; noABHitInd3++){
                  printf(" hit %i - time: %f, evt time: %f offset: %f, energy: %f, core: %u\n",noABHitInd3,noABHitTime(&sortedEvt,noABHitInd3),sortedEvt.header.evtTimeNs,sortedEvt.noABHit[noABHitInd3].timeOffsetNs,sortedEvt.noABHit[noABHitInd3].energy,sortedEvt.noABHit[noABHitInd3].core & 63U);
                  if(noABHitTime(&sortedEvt,noABHitInd3) > 1.0E19){
                    getc(stdin);
                  }
                }
              }
            }*/
            
            hpgeT_hpgeT_le->Fill(sortedEvt.noABHit[noABHitInd2].tsDiff - sortedEvt.noABHit[noABHitInd].tsDiff);
            hpgeT_hpgeT->Fill(tDiff);
            hpgeE_hpgeE->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
            hpgeE_hpgeE->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
            hpgePos_hpgePos->Fill(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U);
            if(fabs(sortedEvt.noABHit[noABHitInd2].energy - sortedEvt.noABHit[noABHitInd].energy) < 0.1){
              hpgePos_hpgePos_lowEDiff->Fill(sortedEvt.noABHit[noABHitInd].core & 63U,sortedEvt.noABHit[noABHitInd2].core & 63U);
            }
            hpge_hpge_dist->Fill(getGeHitDistance(sortedEvt.noABHit[noABHitInd].core & 63U,0,sortedEvt.noABHit[noABHitInd2].core & 63U,0,1)); //FORWARD POSITION (11 cm)
            hpge_hpge_angle->Fill(getGeVector(sortedEvt.noABHit[noABHitInd].core & 63U,0,1).Angle(getGeVector(sortedEvt.noABHit[noABHitInd2].core & 63U,0,1))*180.0/PI); //FORWARD POSITION (11 cm)
            hpgeT_hpgeT_EDiff->Fill(tDiff,sortedEvt.noABHit[noABHitInd2].energy - sortedEvt.noABHit[noABHitInd].energy);
            hpgeT_hpgeT_hpgeE->Fill(tDiff,sortedEvt.noABHit[noABHitInd].energy);
            hpgeT_hpgeT_hpgeE->Fill(tDiff,sortedEvt.noABHit[noABHitInd2].energy);
            if((sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)) && (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6))){
              hpgeT_hpgeT_hpgeE_2CFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd].energy);
              hpgeT_hpgeT_hpgeE_2CFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd2].energy);
            }else if((sortedEvt.noABHit[noABHitInd].core & ((uint8_t)1 << 6)) || (sortedEvt.noABHit[noABHitInd2].core & ((uint8_t)1 << 6))){
              hpgeT_hpgeT_hpgeE_1CFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd].energy);
              hpgeT_hpgeT_hpgeE_1CFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd2].energy);
            }else{
              hpgeT_hpgeT_hpgeE_NoCFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd].energy);
              hpgeT_hpgeT_hpgeE_NoCFDfail->Fill(tDiff,sortedEvt.noABHit[noABHitInd2].energy);
            }
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
            if(getGeVector(sortedEvt.noABHit[noABHitInd].core & 63U,0,1).Angle(getGeVector(sortedEvt.noABHit[noABHitInd2].core & 63U,0,1))*180.0/PI > 175.0){
              if((tDiff >= hpgehpgeTGate[0])&&(tDiff <= hpgehpgeTGate[1])){
                hpgeE_hpgeE_180deg->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
                hpgeE_hpgeE_180deg->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
                if(gate1477 != 0){
                  if((gate1477HitInd != noABHitInd) && (gate1477HitInd != noABHitInd2)){
                    Double_t tDiffCoinc = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,gate1477HitInd);
                    Double_t tDiffCoinc2 = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,gate1477HitInd);
                    if(((tDiffCoinc >= hpgehpgeTGate[0])&&(tDiffCoinc <= hpgehpgeTGate[1])) || ((tDiffCoinc2 >= hpgehpgeTGate[0])&&(tDiffCoinc2 <= hpgehpgeTGate[1]))){
                      hpgeE_hpgeE_180deg_1477gate->Fill(sortedEvt.noABHit[noABHitInd].energy,sortedEvt.noABHit[noABHitInd2].energy);
                      hpgeE_hpgeE_180deg_1477gate->Fill(sortedEvt.noABHit[noABHitInd2].energy,sortedEvt.noABHit[noABHitInd].energy); //symmetrized
                      hpgeE_hpgeE_180deg_sum_1477gate->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy);
                    }
                  }
                }
                hpgeE_hpgeE_180deg_proj->Fill(sortedEvt.noABHit[noABHitInd].energy);
                hpgeE_hpgeE_180deg_proj->Fill(sortedEvt.noABHit[noABHitInd2].energy);
                hpgeE_hpgeE_180deg_sum->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy);
              }
              hpgeE_hpgeE_180deg_sum_tDiff->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy,tDiff);
              if(gate1477 != 0){
                if((gate1477HitInd != noABHitInd) && (gate1477HitInd != noABHitInd2)){
                  Double_t tDiffCoinc = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,gate1477HitInd);
                  Double_t tDiffCoinc2 = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,gate1477HitInd);
                  if(((tDiffCoinc >= hpgehpgeTGate[0])&&(tDiffCoinc <= hpgehpgeTGate[1])) || ((tDiffCoinc2 >= hpgehpgeTGate[0])&&(tDiffCoinc2 <= hpgehpgeTGate[1]))){
                    hpgeE_hpgeE_180deg_sum_tDiff_1477gate->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy,tDiff);
                  }
                }
              }
              if(gate685 != 0){
                if((gate685HitInd != noABHitInd) && (gate685HitInd != noABHitInd2)){
                  Double_t tDiffCoinc = noABHitTime(&sortedEvt,noABHitInd) - noABHitTime(&sortedEvt,gate685HitInd);
                  Double_t tDiffCoinc2 = noABHitTime(&sortedEvt,noABHitInd2) - noABHitTime(&sortedEvt,gate685HitInd);
                  if(((tDiffCoinc >= hpgehpgeTGate[0])&&(tDiffCoinc <= hpgehpgeTGate[1])) || ((tDiffCoinc2 >= hpgehpgeTGate[0])&&(tDiffCoinc2 <= hpgehpgeTGate[1]))){
                    hpgeE_hpgeE_180deg_sum_tDiff_685gate->Fill(sortedEvt.noABHit[noABHitInd].energy + sortedEvt.noABHit[noABHitInd2].energy,tDiff);
                  }
                }
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

    if (jentry % 9713 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << sentries << ", " << 100 * jentry / sentries << "% complete" << "\r" << flush;
  } // analysis tree

  cout << "Entry " << sentries << " of " << sentries << ", 100% complete" << endl;
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
