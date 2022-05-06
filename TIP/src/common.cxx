//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "common.h"

using namespace std;

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
      for(int i=1;i<tip->GetMultiplicity();i++){
        double_t fitTDiff = getTipFitTime(tip->GetTipHit(0),tip_waveform_pretrigger) - getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger);
        if(gate1D(fitTDiff,tiptipTGate[0],tiptipTGate[1])){
          goodTipTipTime = true;
          pass |= (1U<<0);
          pass |= (1U<<i);
        }
      }

      if(goodTipTipTime){
        if(tigress->GetAddbackMultiplicity()>=2){
          //TIGRESS-TIGRESS addback timing
          add_hit = tigress->GetAddbackHit(0);
          suppAdd = add_hit->BGOFired();
          //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
          if(!suppAdd && add_hit->GetEnergy() > 15){
            for(int i=1;i<tigress->GetAddbackMultiplicity();i++){
              add_hit2 = tigress->GetAddbackHit(i);
              suppAdd = add_hit2->BGOFired();
              if(!suppAdd && add_hit2->GetEnergy() > 15){
                double tigtigTDiff = add_hit->GetTime() - add_hit2->GetTime();
                if(gate1D(tigtigTDiff,tigtigTGate[0],tigtigTGate[1])){
                  goodTigTigTime = true;
                  pass |= (1U<<16);
                  pass |= (1U<<(16+i));
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
