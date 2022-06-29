//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "common.h"

using namespace std;

PIDGates::PIDGates(){
  
  //Setup PID gates
  for(int i=0; i<NTIPRING; i++){
    alphaRingCut[i] = new TCutG(Form("ring %i alpha cut",i),8);
    protonRingCut[i] = new TCutG(Form("ring %i proton cut",i),8);
    
    cout << "Creating PID gate for ring " << i << endl;

    switch(i){
      case 9:
        //ring 9
        protonRingCut[i]->SetPoint(0,91,151);
        protonRingCut[i]->SetPoint(1,152,194);
        protonRingCut[i]->SetPoint(2,1460,206);
        protonRingCut[i]->SetPoint(3,1527,172);
        protonRingCut[i]->SetPoint(4,651,155);
        protonRingCut[i]->SetPoint(5,536,149);
        protonRingCut[i]->SetPoint(6,354,149);
        protonRingCut[i]->SetPoint(7,91,151);

        alphaRingCut[i]->SetPoint(0,186,144);
        alphaRingCut[i]->SetPoint(1,543,144);
        alphaRingCut[i]->SetPoint(2,2397,157);
        alphaRingCut[i]->SetPoint(3,2397,133);
        alphaRingCut[i]->SetPoint(4,624,130);
        alphaRingCut[i]->SetPoint(5,368,129);
        alphaRingCut[i]->SetPoint(6,145,131);
        alphaRingCut[i]->SetPoint(7,186,144);
        break;
      case 8:
        //ring 8
        protonRingCut[i]->SetPoint(0,154,137);
        protonRingCut[i]->SetPoint(1,147,190);
        protonRingCut[i]->SetPoint(2,1216,214);
        protonRingCut[i]->SetPoint(3,1497,166);
        protonRingCut[i]->SetPoint(4,730,153);
        protonRingCut[i]->SetPoint(5,552,139);
        protonRingCut[i]->SetPoint(6,319,135);
        protonRingCut[i]->SetPoint(7,154,137);

        alphaRingCut[i]->SetPoint(0,168,135);
        alphaRingCut[i]->SetPoint(1,449,139);
        alphaRingCut[i]->SetPoint(2,2154,143);
        alphaRingCut[i]->SetPoint(3,2154,120);
        alphaRingCut[i]->SetPoint(4,620,115);
        alphaRingCut[i]->SetPoint(5,271,109);
        alphaRingCut[i]->SetPoint(6,100,114);
        alphaRingCut[i]->SetPoint(7,168,135);
        break;
      case 7:
        //ring 7
        protonRingCut[i]->SetPoint(0,160,136);
        protonRingCut[i]->SetPoint(1,229,196);
        protonRingCut[i]->SetPoint(2,2167,215);
        protonRingCut[i]->SetPoint(3,2365,176);
        protonRingCut[i]->SetPoint(4,771,155);
        protonRingCut[i]->SetPoint(5,513,141);
        protonRingCut[i]->SetPoint(6,367,136);
        protonRingCut[i]->SetPoint(7,160,136);

        alphaRingCut[i]->SetPoint(0,237,132);
        alphaRingCut[i]->SetPoint(1,3072,146);
        alphaRingCut[i]->SetPoint(2,3098,115);
        alphaRingCut[i]->SetPoint(3,806,100);
        alphaRingCut[i]->SetPoint(4,367,102);
        alphaRingCut[i]->SetPoint(5,177,106);
        alphaRingCut[i]->SetPoint(6,177,113);
        alphaRingCut[i]->SetPoint(7,237,132);
        break;
      case 6:
        //ring 6
        protonRingCut[i]->SetPoint(0,222,144);
        protonRingCut[i]->SetPoint(1,294,197);
        protonRingCut[i]->SetPoint(2,2385,210);
        protonRingCut[i]->SetPoint(3,2890,178);
        protonRingCut[i]->SetPoint(4,1499,164);
        protonRingCut[i]->SetPoint(5,1077,146);
        protonRingCut[i]->SetPoint(6,459,145);
        protonRingCut[i]->SetPoint(7,222,144);

        alphaRingCut[i]->SetPoint(0,212,142);
        alphaRingCut[i]->SetPoint(1,799,144);
        alphaRingCut[i]->SetPoint(2,3786,152);
        alphaRingCut[i]->SetPoint(3,3838,131);
        alphaRingCut[i]->SetPoint(4,1386,123);
        alphaRingCut[i]->SetPoint(5,510,123);
        alphaRingCut[i]->SetPoint(6,201,127);
        alphaRingCut[i]->SetPoint(7,212,142);
        break;
      case 5:
        //ring 5
        protonRingCut[i]->SetPoint(0,171,136);
        protonRingCut[i]->SetPoint(1,279,192);
        protonRingCut[i]->SetPoint(2,1860,212);
        protonRingCut[i]->SetPoint(3,1918,181);
        protonRingCut[i]->SetPoint(4,896,165);
        protonRingCut[i]->SetPoint(5,664,142);
        protonRingCut[i]->SetPoint(6,432,140);
        protonRingCut[i]->SetPoint(7,171,136);

        alphaRingCut[i]->SetPoint(0,163,134);
        alphaRingCut[i]->SetPoint(1,656,138);
        alphaRingCut[i]->SetPoint(2,2310,151);
        alphaRingCut[i]->SetPoint(3,2302,112);
        alphaRingCut[i]->SetPoint(4,635,103);
        alphaRingCut[i]->SetPoint(5,345,106);
        alphaRingCut[i]->SetPoint(6,113,111);
        alphaRingCut[i]->SetPoint(7,163,134);
        break;
      case 4:
        //ring 4
        protonRingCut[i]->SetPoint(0,189,142);
        protonRingCut[i]->SetPoint(1,329,201);
        protonRingCut[i]->SetPoint(2,1960,215);
        protonRingCut[i]->SetPoint(3,1896,173);
        protonRingCut[i]->SetPoint(4,873,162);
        protonRingCut[i]->SetPoint(5,675,147);
        protonRingCut[i]->SetPoint(6,432,143);
        protonRingCut[i]->SetPoint(7,189,142);

        alphaRingCut[i]->SetPoint(0,176,133);
        alphaRingCut[i]->SetPoint(1,931,141);
        alphaRingCut[i]->SetPoint(2,2670,152);
        alphaRingCut[i]->SetPoint(3,2670,125);
        alphaRingCut[i]->SetPoint(4,675,113);
        alphaRingCut[i]->SetPoint(5,406,105);
        alphaRingCut[i]->SetPoint(6,138,105);
        alphaRingCut[i]->SetPoint(7,176,133);
        break;
      case 3:
        //ring 3
        protonRingCut[i]->SetPoint(0,197,145);
        protonRingCut[i]->SetPoint(1,352,190);
        protonRingCut[i]->SetPoint(2,1125,207);
        protonRingCut[i]->SetPoint(3,2216,210);
        protonRingCut[i]->SetPoint(4,2332,181);
        protonRingCut[i]->SetPoint(5,963,161);
        protonRingCut[i]->SetPoint(6,630,146);
        protonRingCut[i]->SetPoint(7,197,145);

        alphaRingCut[i]->SetPoint(0,213,134);
        alphaRingCut[i]->SetPoint(1,1195,142);
        alphaRingCut[i]->SetPoint(2,3198,150);
        alphaRingCut[i]->SetPoint(3,3136,126);
        alphaRingCut[i]->SetPoint(4,824,116);
        alphaRingCut[i]->SetPoint(5,499,106);
        alphaRingCut[i]->SetPoint(6,159,109);
        alphaRingCut[i]->SetPoint(7,213,134);
        break;
      case 2:
        //ring 2
        protonRingCut[i]->SetPoint(0,204,150);
        protonRingCut[i]->SetPoint(1,304,186);
        protonRingCut[i]->SetPoint(2,2118,221);
        protonRingCut[i]->SetPoint(3,2276,184);
        protonRingCut[i]->SetPoint(4,1119,165);
        protonRingCut[i]->SetPoint(5,695,149);
        protonRingCut[i]->SetPoint(6,454,146);
        protonRingCut[i]->SetPoint(7,204,150);

        alphaRingCut[i]->SetPoint(0,237,138);
        alphaRingCut[i]->SetPoint(1,1319,144);
        alphaRingCut[i]->SetPoint(2,3491,153);
        alphaRingCut[i]->SetPoint(3,3516,131);
        alphaRingCut[i]->SetPoint(4,903,111);
        alphaRingCut[i]->SetPoint(5,529,102);
        alphaRingCut[i]->SetPoint(6,188,104);
        alphaRingCut[i]->SetPoint(7,237,138);
        break;
      case 1:
        //ring 1
        protonRingCut[i]->SetPoint(0,810,176);
        protonRingCut[i]->SetPoint(1,1073,212);
        protonRingCut[i]->SetPoint(2,3610,216);
        protonRingCut[i]->SetPoint(3,4004,191);
        protonRingCut[i]->SetPoint(4,2618,155);
        protonRingCut[i]->SetPoint(5,1554,149);
        protonRingCut[i]->SetPoint(6,1102,149);
        protonRingCut[i]->SetPoint(7,810,176);

        alphaRingCut[i]->SetPoint(0,1044,138);
        alphaRingCut[i]->SetPoint(1,2662,141);
        alphaRingCut[i]->SetPoint(2,4849,149);
        alphaRingCut[i]->SetPoint(3,4631,129);
        alphaRingCut[i]->SetPoint(4,1904,119);
        alphaRingCut[i]->SetPoint(5,1466,112);
        alphaRingCut[i]->SetPoint(6,971,114);
        alphaRingCut[i]->SetPoint(7,1044,138);
        break;
      case 0:
      default:
        //ring 0
        protonRingCut[i]->SetPoint(0,493,150);
        protonRingCut[i]->SetPoint(1,622,227);
        protonRingCut[i]->SetPoint(2,2243,236);
        protonRingCut[i]->SetPoint(3,3048,238);
        protonRingCut[i]->SetPoint(4,2992,173);
        protonRingCut[i]->SetPoint(5,1572,157);
        protonRingCut[i]->SetPoint(6,1192,150);
        protonRingCut[i]->SetPoint(7,493,150);

        alphaRingCut[i]->SetPoint(0,532,136);
        alphaRingCut[i]->SetPoint(1,1057,142);
        alphaRingCut[i]->SetPoint(2,2276,152);
        alphaRingCut[i]->SetPoint(3,3820,155);
        alphaRingCut[i]->SetPoint(4,3719,123);
        alphaRingCut[i]->SetPoint(5,1516,117);
        alphaRingCut[i]->SetPoint(6,509,104);
        alphaRingCut[i]->SetPoint(7,532,136);
        break;
    }

  }
  

  
  
}

//get the particle type for a tip hit
//0=unidentified
//1=proton
//4=alpha
Int_t getParticleType(TTipHit *tip_hit, PIDGates *gates){
  double_t tipPID = tip_hit->GetPID();

  //check if the hit is a proton or alpha
  if(tipPID>=0){ //PID was found
    Int_t ring = getTIPRing(tip_hit->GetTipChannel());
    if(gates->protonRingCut[ring]->IsInside(tip_hit->GetEnergy(),tipPID)){
      return 1;
    }else if(gates->alphaRingCut[ring]->IsInside(tip_hit->GetEnergy(),tipPID)){
      return 4;
    }
  }
  return 0;
}

//given a TIGRESS hit and the (calibrated in MeV) TIP hits in a fusion-evaporation
//event, determine the Doppler corrected TIGRESS energy
double_t getEDoppFusEvap(TTigressHit *add_hit, TTip *tip, const uint64_t passedtimeGate, PIDGates *gates){
  double_t resM = compoundM_AMU*AMU; //residual mass prior to particle evaporation
  TVector3 p_compound(0,0,resM*betaCompound/(1.0-betaCompound*betaCompound));
  //cout << "p_compound: " << p_compound.X() << " " << p_compound.Y() << " " << p_compound.Z() << endl;
  TVector3 p_part;
  TTipHit *tip_hit;
  double_t resBeta = 0.;
  if(tip->GetMultiplicity()<=MAXNUMTIPHIT){
    for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
      if(passedtimeGate&(1ULL<<tipHitInd)){
        tip_hit = tip->GetTipHit(tipHitInd);
        switch(getParticleType(tip_hit,gates)){ //see common.cxx
          case 4:
            p_part = tip_hit->GetPosition();
            //cout << "p_alpha " << tipHitInd << ": " << p_part.X() << " " << p_part.Y() << " " << p_part.Z() << endl;
            //p_part.SetMag(sqrt(2.0*4.0015*AMU*tip_hit->GetEnergy()));
            p_part.SetMag(sqrt(2.0*4.0015*AMU*10)); //TEMPORARY - ASSUME 10 MeV
            p_compound -= p_part;
            resM -= 2.0*4.0015*AMU; //technically not true due to binding energy!
            break;
          case 1:
            p_part = tip_hit->GetPosition();
            //cout << "p_proton " << tipHitInd << ": " << p_part.X() << " " << p_part.Y() << " " << p_part.Z() << endl;
            //p_part.SetMag(sqrt(2.0*938.272*tip_hit->GetEnergy()));
            p_part.SetMag(sqrt(2.0*938.272*10)); //TEMPORARY - ASSUME 10 MeV
            p_compound -= p_part;
            resM -= 2.0*938.272; //technically not true due to binding energy!
            break;
          case 0:
          default:
            break;
        }
      }
    }
  }
  //cout << "p_res: " << p_compound.X() << " " << p_compound.Y() << " " << p_compound.Z() << endl;
  resBeta = sqrt(1.0 - 1.0/(1.0 + pow(p_compound.Mag()/resM,2)));
  //cout << "resBeta: " << resBeta << ", E: " << add_hit->GetDoppler(resBeta,&p_compound) << endl;
  return add_hit->GetDoppler(resBeta,&p_compound);
  
}

double_t getTipFitTime(TTipHit *tip_hit, const Int_t pretrigger_samples){
  
  double_t tTipFit = 0.;
  if((tip_hit->GetPID() > -1000.)&&(tip_hit->GetWaveform()->size() > 50)){
    tTipFit = tip_hit->GetFitTime() * 10.; //fit time in samples
    tTipFit += ((tip_hit->GetTimeStamp()) + gRandom->Uniform()) * tip_hit->GetTimeStampUnit();
    tTipFit -= pretrigger_samples * 10.;
  }
  //cout << "fit time: " << tTipFit << ", CFD time: " << tip_hit->GetTime() << endl;
  return tTipFit;
}

//function defines the Compton suppression timing window for TIGRESS
bool ExptSuppression(TDetectorHit* tig, TBgoHit& bgo){
   Int_t dCfd = static_cast<TTigressHit*>(tig)->GetCfd() - bgo.GetCfd();
   return ((dCfd > tigBGOTGate[0] && dCfd < -tigBGOTGate[1]) && (tig->GetDetector() == bgo.GetDetector()) && (bgo.GetEnergy() > 0));
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
//bits 0 to (MAXNUMTIPHIT-1): TIP hit indices
//bits MAXNUMTIPHIT to MAXNUMTIPHIT+MAXNUMTIGHIT: TIGRESS hit indices
uint64_t passesTimeGate(TTigress *tigress, TTip *tip){

  //Defining Pointers
  TTigressHit *add_hit, *add_hit2;
  bool suppAdd = false;

  bool goodTipTipTime = false;
  bool goodTigTigTime = false;
  bool goodTipTigTime = false;
  uint64_t pass = 0;

  //evaluate timing conditions
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
          if(tip->GetTipHit(j)->GetKValue() != noPileupKValue){
            //pileup
            continue;
          }
          double_t fitT1 = getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger);
          double_t fitT2 = getTipFitTime(tip->GetTipHit(j),tip_waveform_pretrigger);
          if((fitT1 != 0.)&&(fitT2 != 0.)){
            double_t fitTDiff = fitT1 - fitT2;
            if(gate1D(fitTDiff,tiptipTGate[i],tiptipTGate[j])){
              goodTipTipTime = true;
              pass |= (1ULL<<i);
              pass |= (1ULL<<j);
            }
          }
        }
      }

      if(goodTipTipTime){
        pass |= (1ULL<<61);
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
                if(add_hit2->GetKValue() != noPileupKValue){
                  //pileup
                  continue;
                }
                suppAdd = add_hit2->BGOFired();
                if(!suppAdd && add_hit2->GetEnergy() > 15){
                  double tigtigTDiff = add_hit->GetTime() - add_hit2->GetTime();
                  if(gate1D(tigtigTDiff,tigtigTGate[0],tigtigTGate[1])){
                    goodTigTigTime = true;
                    pass |= (1ULL<<(MAXNUMTIPHIT+i));
                    pass |= (1ULL<<(MAXNUMTIPHIT+j));
                  }
                }
              }
            }
          }
        }
      }

      if(goodTigTigTime){
        pass |= (1ULL<<62);
        for(int tipHitInd=0;tipHitInd<tip->GetMultiplicity();tipHitInd++){
          if(pass&(1ULL<<tipHitInd)){
            for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
              if(pass&(1ULL<<(MAXNUMTIPHIT+tigHitIndAB))){
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

  if(goodTipTigTime){
    pass |= (1ULL<<63);
    return pass;
  }else{
    pass = 0;
    if(goodTipTipTime){
      pass |= (1ULL<<61);
    }
    if(goodTigTigTime){
      pass |= (1ULL<<62);
    }
    return pass;
  }
  
}
