//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "common.h"

using namespace std;

void setupPIDGates(){
  //Setup PID gates
  printf("Creating PID gates\n");
  //individual gates for each detector
  for(int i=0; i<NTIPRING; i++){
    alphaRingCut[i] = new TCutG(Form("ring %i alpha cut",i),7);
    protonRingCut[i] = new TCutG(Form("ring %i proton cut",i),7);
  }

  //ring 0
  protonRingCut[0]->SetPoint(0,10.4,148);
  protonRingCut[0]->SetPoint(1,8.8,212);
  protonRingCut[0]->SetPoint(2,17.4,257);
  protonRingCut[0]->SetPoint(3,31,273);
  protonRingCut[0]->SetPoint(4,31,191);
  protonRingCut[0]->SetPoint(5,22.3,173);
  protonRingCut[0]->SetPoint(6,17.3,144);

  alphaRingCut[0]->SetPoint(0,21.3,155);
  alphaRingCut[0]->SetPoint(1,48,176);
  alphaRingCut[0]->SetPoint(2,47,125);
  alphaRingCut[0]->SetPoint(3,29.2,115);
  alphaRingCut[0]->SetPoint(4,23.8,107);
  alphaRingCut[0]->SetPoint(5,20.3,112);
  alphaRingCut[0]->SetPoint(6,18.9,128);

  //ring 1
  protonRingCut[1]->SetPoint(0,10.4,128);
  protonRingCut[1]->SetPoint(1,8.8,202);
  protonRingCut[1]->SetPoint(2,17.4,247);
  protonRingCut[1]->SetPoint(3,31,253);
  protonRingCut[1]->SetPoint(4,31,181);
  protonRingCut[1]->SetPoint(5,22.3,163);
  protonRingCut[1]->SetPoint(6,17.3,124);

  alphaRingCut[1]->SetPoint(0,21.3,150);
  alphaRingCut[1]->SetPoint(1,48,176);
  alphaRingCut[1]->SetPoint(2,47,125);
  alphaRingCut[1]->SetPoint(3,29.2,115);
  alphaRingCut[1]->SetPoint(4,23.8,107);
  alphaRingCut[1]->SetPoint(5,20.3,112);
  alphaRingCut[1]->SetPoint(6,18.9,128);

  //ring 2
  protonRingCut[2]->SetPoint(0,7.2,134);
  protonRingCut[2]->SetPoint(1,9.2,241);
  protonRingCut[2]->SetPoint(2,21.8,262);
  protonRingCut[2]->SetPoint(3,31,237);
  protonRingCut[2]->SetPoint(4,30.5,193);
  protonRingCut[2]->SetPoint(5,18.5,156);
  protonRingCut[2]->SetPoint(6,13,128);

  alphaRingCut[2]->SetPoint(0,15.8,129);
  alphaRingCut[2]->SetPoint(1,19.5,154);
  alphaRingCut[2]->SetPoint(2,54.3,174);
  alphaRingCut[2]->SetPoint(3,53,134);
  alphaRingCut[2]->SetPoint(4,28.3,120);
  alphaRingCut[2]->SetPoint(5,18.6,106);
  alphaRingCut[2]->SetPoint(6,15.4,112);

  //ring 3
  protonRingCut[3]->SetPoint(0,6.9,124);
  protonRingCut[3]->SetPoint(1,7.9,227);
  protonRingCut[3]->SetPoint(2,17.2,247);
  protonRingCut[3]->SetPoint(3,27.9,247);
  protonRingCut[3]->SetPoint(4,29,206.7);
  protonRingCut[3]->SetPoint(5,18,161.3);
  protonRingCut[3]->SetPoint(6,11.5,126.7);

  alphaRingCut[3]->SetPoint(0,13.1,126);
  alphaRingCut[3]->SetPoint(1,14.2,143.7);
  alphaRingCut[3]->SetPoint(2,21.7,156);
  alphaRingCut[3]->SetPoint(3,45.4,161);
  alphaRingCut[3]->SetPoint(4,45.3,134);
  alphaRingCut[3]->SetPoint(5,22,115.2);
  alphaRingCut[3]->SetPoint(6,13,104.5);

  //ring 4
  protonRingCut[4]->SetPoint(0,6.7,119.8);
  protonRingCut[4]->SetPoint(1,7.9,219);
  protonRingCut[4]->SetPoint(2,20.5,250.6);
  protonRingCut[4]->SetPoint(3,28.3,217.5);
  protonRingCut[4]->SetPoint(4,25.6,186);
  protonRingCut[4]->SetPoint(5,13.7,142.9);
  protonRingCut[4]->SetPoint(6,8.4,127.5);

  alphaRingCut[4]->SetPoint(0,12.1,115.2);
  alphaRingCut[4]->SetPoint(1,15.1,140.6);
  alphaRingCut[4]->SetPoint(2,46.5,170.6);
  alphaRingCut[4]->SetPoint(3,47.4,131.4);
  alphaRingCut[4]->SetPoint(4,19.1,109);
  alphaRingCut[4]->SetPoint(5,14.3,97.5);
  alphaRingCut[4]->SetPoint(6,12.4,109.1);

  //ring 5
  protonRingCut[5]->SetPoint(0,6.1,137.5);
  protonRingCut[5]->SetPoint(1,7.4,209.8);
  protonRingCut[5]->SetPoint(2,25,251.4);
  protonRingCut[5]->SetPoint(3,26.3,191);
  protonRingCut[5]->SetPoint(4,16.5,159);
  protonRingCut[5]->SetPoint(5,12.2,138);
  protonRingCut[5]->SetPoint(6,9.7,136);

  alphaRingCut[5]->SetPoint(0,12.2,132.9);
  alphaRingCut[5]->SetPoint(1,17.3,148.3);
  alphaRingCut[5]->SetPoint(2,38.4,150.6);
  alphaRingCut[5]->SetPoint(3,39.6,119.8);
  alphaRingCut[5]->SetPoint(4,17.4,107.5);
  alphaRingCut[5]->SetPoint(5,11.1,103.6);
  alphaRingCut[5]->SetPoint(6,11.3,118.3);

  //ring 6
  protonRingCut[6]->SetPoint(0,0,0);
  protonRingCut[6]->SetPoint(1,0,0);
  protonRingCut[6]->SetPoint(2,0,0);
  protonRingCut[6]->SetPoint(3,0,0);
  protonRingCut[6]->SetPoint(4,0,0);
  protonRingCut[6]->SetPoint(5,0,0);
  protonRingCut[6]->SetPoint(6,0,0);

  alphaRingCut[6]->SetPoint(0,0,0);
  alphaRingCut[6]->SetPoint(1,0,0);
  alphaRingCut[6]->SetPoint(2,0,0);
  alphaRingCut[6]->SetPoint(3,0,0);
  alphaRingCut[6]->SetPoint(4,0,0);
  alphaRingCut[6]->SetPoint(5,0,0);
  alphaRingCut[6]->SetPoint(6,0,0);

  //ring 7
  protonRingCut[7]->SetPoint(0,6.3,158);
  protonRingCut[7]->SetPoint(1,7.1,216);
  protonRingCut[7]->SetPoint(2,22.4,239);
  protonRingCut[7]->SetPoint(3,23.3,180);
  protonRingCut[7]->SetPoint(4,15.7,156);
  protonRingCut[7]->SetPoint(5,13.2,144);
  protonRingCut[7]->SetPoint(6,9.9,142);

  alphaRingCut[7]->SetPoint(0,12.4,102);
  alphaRingCut[7]->SetPoint(1,14.5,143);
  alphaRingCut[7]->SetPoint(2,27.4,156);
  alphaRingCut[7]->SetPoint(3,28,108);
  alphaRingCut[7]->SetPoint(4,14.4,98.3);
  alphaRingCut[7]->SetPoint(5,12.9,107);
  alphaRingCut[7]->SetPoint(6,12.6,105);

  //ring 8
  protonRingCut[8]->SetPoint(0,5.8,156);
  protonRingCut[8]->SetPoint(1,6.7,204);
  protonRingCut[8]->SetPoint(2,17.6,225);
  protonRingCut[8]->SetPoint(3,22,196);
  protonRingCut[8]->SetPoint(4,21.4,171.3);
  protonRingCut[8]->SetPoint(5,14.1,152);
  protonRingCut[8]->SetPoint(6,11.3,141.4);

  alphaRingCut[8]->SetPoint(0,13,131);
  alphaRingCut[8]->SetPoint(1,16.1,139);
  alphaRingCut[8]->SetPoint(2,25.6,137.5);
  alphaRingCut[8]->SetPoint(3,26,103);
  alphaRingCut[8]->SetPoint(4,14.1,96.5);
  alphaRingCut[8]->SetPoint(5,11.5,105);
  alphaRingCut[8]->SetPoint(6,12.1,124);

  //ring 9
  protonRingCut[9]->SetPoint(0,5.8,156);
  protonRingCut[9]->SetPoint(1,6.7,204);
  protonRingCut[9]->SetPoint(2,17.6,225);
  protonRingCut[9]->SetPoint(3,22,196);
  protonRingCut[9]->SetPoint(4,21.4,171.3);
  protonRingCut[9]->SetPoint(5,14.1,152);
  protonRingCut[9]->SetPoint(6,11.3,141.4);

  alphaRingCut[9]->SetPoint(0,13,131);
  alphaRingCut[9]->SetPoint(1,16.1,139);
  alphaRingCut[9]->SetPoint(2,25.6,137.5);
  alphaRingCut[9]->SetPoint(3,26,103);
  alphaRingCut[9]->SetPoint(4,14.1,96.5);
  alphaRingCut[9]->SetPoint(5,11.5,105);
  alphaRingCut[9]->SetPoint(6,12.1,124);
  
}

double_t getTipFitTime(TTipHit *tip_hit, const Int_t pretrigger_samples){
  
  double_t tTipFit = 0.;
  const std::vector<Short_t> *wf = tip_hit->GetWaveform();
  TPulseAnalyzer pulse;
  pulse.SetData(*wf, 0);
  if(pulse.CsIPID() > -1000.0){
    if(wf->size() > 50){
      tTipFit = pulse.CsIt0() * 10.;
      tTipFit += ((tip_hit->GetTimeStamp()) + gRandom->Uniform()) * tip_hit->GetTimeStampUnit();
      tTipFit -= pretrigger_samples * 10.;
    }
  }
  
  //cout << "fit time: " << tTipFit << endl;
  return tTipFit;
}

//function defines the Compton suppression timing window for TIGRESS
bool S1232Suppression(TDetectorHit* tig, TBgoHit& bgo){
   Int_t dCfd = static_cast<TTigressHit*>(tig)->GetCfd() - bgo.GetCfd();
   return ((dCfd > -375 && dCfd < -100) && (tig->GetDetector() == bgo.GetDetector()) && (bgo.GetEnergy() > 0));
}


bool gate1D(const Double_t value, const Double_t min, const Double_t max){
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


//returns a bit-pattern specifying which TIP and TIGRESS events
//passed the timing gates specified in common.h
//bits 0-15: TIP hit indices
//bits 16-31: TIGRESS hit indices
uint32_t passesTimeGate(TTigress *tigress, TTip *tip){

  //Defining Pointers
  TTigressHit *add_hit, *add_hit2;
  bool suppAdd = false;

  bool goodTipTipTime = false;
  bool goodTigTigTime = false;
  bool goodTipTigTime = false;
  uint32_t pass = 0;

  //first, evaluate timing conditions
  if(tip && tigress){

    if(tip->GetMultiplicity()>MAXNUMTIPHIT){
      return 0;
    }
    if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT){
      return 0;
    }
    
    if(tip->GetMultiplicity()>=2){
      
      //TIP-TIP timing
      for(int i=0;i<tip->GetMultiplicity();i++){
        if(tip->GetTipHit(i)->GetKValue() != noPileupKValue){
          //pileup
          continue;
        }
        for(int j=i+1;j<tip->GetMultiplicity();j++){
          double_t fitTDiff = getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger) - getTipFitTime(tip->GetTipHit(j),tip_waveform_pretrigger);
          if(gate1D(fitTDiff,tiptipTGate[i],tiptipTGate[j])){
            goodTipTipTime = true;
            pass |= (1U<<i);
            pass |= (1U<<j);
          }
        }
      }

      if(goodTipTipTime){
        if(tigress->GetAddbackMultiplicity()>=2){
          for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
            //TIGRESS-TIGRESS addback timing
            add_hit = tigress->GetAddbackHit(i);
            if(add_hit->GetKValue() != noPileupKValue){
              //pileup
              continue;
            }
            suppAdd = add_hit->BGOFired();
            //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
            if(!suppAdd && add_hit->GetEnergy() > 15){
              for(int j=i+1;j<tigress->GetAddbackMultiplicity();j++){
                add_hit2 = tigress->GetAddbackHit(j);
                suppAdd = add_hit2->BGOFired();
                if(!suppAdd && add_hit2->GetEnergy() > 15){
                  double tigtigTDiff = add_hit->GetTime() - add_hit2->GetTime();
                  if(gate1D(tigtigTDiff,tigtigTGate[0],tigtigTGate[1])){
                    goodTigTigTime = true;
                    pass |= (1U<<(16+i));
                    pass |= (1U<<(16+j));
                  }
                }
              }
            }
          }
        }
      }

      if(goodTigTigTime){
        for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
          if(pass&(1U<<tipHitInd)){
            for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
              if(pass&(1U<<(tigHitIndAB+16))){
                add_hit = tigress->GetAddbackHit(tigHitIndAB);
                suppAdd = add_hit->BGOFired();
                if(!suppAdd && add_hit->GetEnergy() > 15){
                  double_t tiptigTDiff = getTipFitTime(tip->GetTipHit(tipHitInd),tip_waveform_pretrigger) - add_hit->GetTime();
                  if(gate1D(tiptigTDiff,tiptigTGate[0],tiptigTGate[1])){
                    goodTipTigTime = true;
                    break;
                  }
                }
              }
            }
            if(goodTipTigTime){
              break;
            }
          }
        }
      }

    }

  }

  if((!goodTigTigTime)||(!goodTipTipTime)||(!goodTipTigTime)){
    return 0;
  }

  return pass;
}
