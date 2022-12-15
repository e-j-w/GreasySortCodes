//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "../include/common.h"

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
        protonRingCut[i]->SetPoint(0,2,51);
        protonRingCut[i]->SetPoint(1,2.5,94);
        protonRingCut[i]->SetPoint(2,13.5,104.5);
        protonRingCut[i]->SetPoint(3,13.3,82.9);
        protonRingCut[i]->SetPoint(4,7.4,50.1);
        protonRingCut[i]->SetPoint(5,5.6,29.3);
        protonRingCut[i]->SetPoint(6,3.7,27.7);
        protonRingCut[i]->SetPoint(7,2,51);

        alphaRingCut[i]->SetPoint(0,2.2,24);
        alphaRingCut[i]->SetPoint(1,5,24);
        alphaRingCut[i]->SetPoint(2,20,37);
        alphaRingCut[i]->SetPoint(3,20,13);
        alphaRingCut[i]->SetPoint(4,5.5,8);
        alphaRingCut[i]->SetPoint(5,3.7,2.5);
        alphaRingCut[i]->SetPoint(6,2,11);
        alphaRingCut[i]->SetPoint(7,2.2,24);
        break;
      case 8:
        //ring 8
        protonRingCut[i]->SetPoint(0,1,37);
        protonRingCut[i]->SetPoint(1,3.5,94.9);
        protonRingCut[i]->SetPoint(2,12.4,112);
        protonRingCut[i]->SetPoint(3,12.7,83.6);
        protonRingCut[i]->SetPoint(4,7.1,47.8);
        protonRingCut[i]->SetPoint(5,5.7,28.5);
        protonRingCut[i]->SetPoint(6,3.6,24.1);
        protonRingCut[i]->SetPoint(7,1,37);

        alphaRingCut[i]->SetPoint(0,1.2,22.8);
        alphaRingCut[i]->SetPoint(1,4.9,23.7);
        alphaRingCut[i]->SetPoint(2,19,27);
        alphaRingCut[i]->SetPoint(3,19,20);
        alphaRingCut[i]->SetPoint(4,6.3,8.9);
        alphaRingCut[i]->SetPoint(5,4.2,3);
        alphaRingCut[i]->SetPoint(6,2.7,6.5);
        alphaRingCut[i]->SetPoint(7,1.2,22.8);
        break;
      case 7:
        //ring 7
        protonRingCut[i]->SetPoint(0,2.5,33.5);
        protonRingCut[i]->SetPoint(1,2,96);
        protonRingCut[i]->SetPoint(2,19,115);
        protonRingCut[i]->SetPoint(3,20,76);
        protonRingCut[i]->SetPoint(4,9.8,60);
        protonRingCut[i]->SetPoint(5,7.6,45);
        protonRingCut[i]->SetPoint(6,6.75,37);
        protonRingCut[i]->SetPoint(7,2.5,33.5);

        alphaRingCut[i]->SetPoint(0,1.9,26);
        alphaRingCut[i]->SetPoint(1,26,46);
        alphaRingCut[i]->SetPoint(2,27,15);
        alphaRingCut[i]->SetPoint(3,7,0);
        alphaRingCut[i]->SetPoint(4,3,2);
        alphaRingCut[i]->SetPoint(5,1.5,6);
        alphaRingCut[i]->SetPoint(6,1.2,13);
        alphaRingCut[i]->SetPoint(7,1.9,26);
        break;
      case 6:
        //ring 6
        protonRingCut[i]->SetPoint(0,1.8,31);
        protonRingCut[i]->SetPoint(1,2.5,97);
        protonRingCut[i]->SetPoint(2,20,110);
        protonRingCut[i]->SetPoint(3,24,78);
        protonRingCut[i]->SetPoint(4,12,64);
        protonRingCut[i]->SetPoint(5,9,46);
        protonRingCut[i]->SetPoint(6,7,37.2);
        protonRingCut[i]->SetPoint(7,1.8,31);

        alphaRingCut[i]->SetPoint(0,1.8,22);
        alphaRingCut[i]->SetPoint(1,7,24);
        alphaRingCut[i]->SetPoint(2,32,32);
        alphaRingCut[i]->SetPoint(3,33,11);
        alphaRingCut[i]->SetPoint(4,10,3);
        alphaRingCut[i]->SetPoint(5,4.2,3);
        alphaRingCut[i]->SetPoint(6,1.5,7);
        alphaRingCut[i]->SetPoint(7,1.8,22);
        break;
      case 5:
        //ring 5
        protonRingCut[i]->SetPoint(0,1.2,36);
        protonRingCut[i]->SetPoint(1,2.5,92);
        protonRingCut[i]->SetPoint(2,19.5,118);
        protonRingCut[i]->SetPoint(3,20,90);
        protonRingCut[i]->SetPoint(4,11.4,69.3);
        protonRingCut[i]->SetPoint(5,8.7,51);
        protonRingCut[i]->SetPoint(6,8.3,44.8);
        protonRingCut[i]->SetPoint(7,1.2,36);

        alphaRingCut[i]->SetPoint(0,1.2,34);
        alphaRingCut[i]->SetPoint(1,5.8,38);
        alphaRingCut[i]->SetPoint(2,20,51);
        alphaRingCut[i]->SetPoint(3,20,12);
        alphaRingCut[i]->SetPoint(4,5.5,3);
        alphaRingCut[i]->SetPoint(5,2.6,3.1);
        alphaRingCut[i]->SetPoint(6,1,11);
        alphaRingCut[i]->SetPoint(7,1.2,34);
        break;
      case 4:
        //ring 4
        protonRingCut[i]->SetPoint(0,1.5,42);
        protonRingCut[i]->SetPoint(1,2.8,101);
        protonRingCut[i]->SetPoint(2,16,115);
        protonRingCut[i]->SetPoint(3,14.3,91.4);
        protonRingCut[i]->SetPoint(4,11.8,71.2);
        protonRingCut[i]->SetPoint(5,10,58.3);
        protonRingCut[i]->SetPoint(6,8.9,48.2);
        protonRingCut[i]->SetPoint(7,1.5,42);

        alphaRingCut[i]->SetPoint(0,1.2,33);
        alphaRingCut[i]->SetPoint(1,8,41);
        alphaRingCut[i]->SetPoint(2,23,52);
        alphaRingCut[i]->SetPoint(3,23,25);
        alphaRingCut[i]->SetPoint(4,10.6,14);
        alphaRingCut[i]->SetPoint(5,7.27,3);
        alphaRingCut[i]->SetPoint(6,1,5);
        alphaRingCut[i]->SetPoint(7,1.2,33);
        break;
      case 3:
        //ring 3
        protonRingCut[i]->SetPoint(0,1.5,45);
        protonRingCut[i]->SetPoint(1,3,90);
        protonRingCut[i]->SetPoint(2,9.6,112);
        protonRingCut[i]->SetPoint(3,18.1,120);
        protonRingCut[i]->SetPoint(4,19.6,94.8);
        protonRingCut[i]->SetPoint(5,13.6,79.7);
        protonRingCut[i]->SetPoint(6,11.6,50);
        protonRingCut[i]->SetPoint(7,1.5,45);

        alphaRingCut[i]->SetPoint(0,2.9,37.6);
        alphaRingCut[i]->SetPoint(1,10,42);
        alphaRingCut[i]->SetPoint(2,28,50);
        alphaRingCut[i]->SetPoint(3,28,26);
        alphaRingCut[i]->SetPoint(4,12.4,14.3);
        alphaRingCut[i]->SetPoint(5,9.6,5.1);
        alphaRingCut[i]->SetPoint(6,1.5,2.9);
        alphaRingCut[i]->SetPoint(7,2.9,37.6);
        break;
      case 2:
        //ring 2
        protonRingCut[i]->SetPoint(0,1.5,50);
        protonRingCut[i]->SetPoint(1,2.6,95);
        protonRingCut[i]->SetPoint(2,19.3,124.1);
        protonRingCut[i]->SetPoint(3,23,97.7);
        protonRingCut[i]->SetPoint(4,16,71);
        protonRingCut[i]->SetPoint(5,11,54);
        protonRingCut[i]->SetPoint(6,7,46);
        protonRingCut[i]->SetPoint(7,1.5,50);

        alphaRingCut[i]->SetPoint(0,1.5,38);
        alphaRingCut[i]->SetPoint(1,11,44);
        alphaRingCut[i]->SetPoint(2,30,53);
        alphaRingCut[i]->SetPoint(3,31,31);
        alphaRingCut[i]->SetPoint(4,12.7,9.7);
        alphaRingCut[i]->SetPoint(5,4,2);
        alphaRingCut[i]->SetPoint(6,1,4);
        alphaRingCut[i]->SetPoint(7,1.5,38);
        break;
      case 1:
        //ring 1
        protonRingCut[i]->SetPoint(0,3,76);
        protonRingCut[i]->SetPoint(1,8,112);
        protonRingCut[i]->SetPoint(2,30,116);
        protonRingCut[i]->SetPoint(3,35,91);
        protonRingCut[i]->SetPoint(4,22,55);
        protonRingCut[i]->SetPoint(5,12,49);
        protonRingCut[i]->SetPoint(6,2.1,49);
        protonRingCut[i]->SetPoint(7,3,76);

        alphaRingCut[i]->SetPoint(0,2.7,38);
        alphaRingCut[i]->SetPoint(1,22,45.5);
        alphaRingCut[i]->SetPoint(2,40,49);
        alphaRingCut[i]->SetPoint(3,41,29);
        alphaRingCut[i]->SetPoint(4,21,13);
        alphaRingCut[i]->SetPoint(5,15.6,3);
        alphaRingCut[i]->SetPoint(6,2.4,3.6);
        alphaRingCut[i]->SetPoint(7,2.7,38);
        break;
      case 0:
      default:
        //ring 0
        protonRingCut[i]->SetPoint(0,4,50);
        protonRingCut[i]->SetPoint(1,5,127);
        protonRingCut[i]->SetPoint(2,19,136);
        protonRingCut[i]->SetPoint(3,26,138);
        protonRingCut[i]->SetPoint(4,26,73);
        protonRingCut[i]->SetPoint(5,12,57);
        protonRingCut[i]->SetPoint(6,9,50);
        protonRingCut[i]->SetPoint(7,4,50);

        alphaRingCut[i]->SetPoint(0,3,27);
        alphaRingCut[i]->SetPoint(1,8,42);
        alphaRingCut[i]->SetPoint(2,18,52);
        alphaRingCut[i]->SetPoint(3,32,55);
        alphaRingCut[i]->SetPoint(4,30,23);
        alphaRingCut[i]->SetPoint(5,16.1,6.7);
        alphaRingCut[i]->SetPoint(6,5,4);
        alphaRingCut[i]->SetPoint(7,3,27);
        break;
    }

  }
  
}

//get the particle type for a tip hit
//0=unidentified
//1=proton
//4=alpha
Int_t getParticleTypePID(double_t tipPID, double_t energy, Int_t detNum, PIDGates *gates){
  //check if the hit is a proton or alpha
  if(tipPID>=0){ //PID was found
    Int_t ring = getTIPRing(detNum);
    if((ring >= 0)&&(ring < NTIPRING)){
      if(gates->protonRingCut[ring]->IsInside(energy,tipPID)){
        return 1;
      }else if(gates->alphaRingCut[ring]->IsInside(energy,tipPID)){
        return 4;
      }
    }
  }
  return 0;
}
Int_t getParticleType(TTipHit *tip_hit, PIDGates *gates){
  double_t tipPID = tip_hit->GetPID();

  //check if the hit is a proton or alpha
  return getParticleTypePID(tipPID,tip_hit->GetEnergy(),tip_hit->GetTipChannel(),gates);
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
  //gets the fit time, without relying on random number generation from GRSISort
  
  double_t tTipFit = 0.;
  //if((tip_hit->GetPID() > -1000.)&&(tip_hit->GetWaveform()->size() > 50)){
  if(tip_hit->GetPID() > 0.){
    tTipFit = tip_hit->GetFitTime() * 10.; //fit time in samples
    tTipFit += tip_hit->GetTimeStamp() * tip_hit->GetTimeStampUnit();
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

Int_t getTIGRESSSegmentRing(const float theta){
  if((theta < 0.0f)||(theta > 180.0f)){
    cout << "Invalid TIGRESS angle: " << theta << endl;
    return 0;
  }else if(theta <= 35.0f){
    return 0;
  }else if(theta <= 45.0f){
    return 1;
  }else if(theta <= 55.0f){
    return 2;
  }else if(theta <= 67.5f){
    return 3;
  }else if(theta <= 80.0f){
    return 4;
  }else if(theta <= 90.0f){
    return 5;
  }else if(theta <= 100.0f){
    return 6;
  }else if(theta <= 112.5f){
    return 7;
  }else if(theta <= 125.0f){
    return 8;
  }else if(theta <= 135.0f){
    return 9;
  }else if(theta <= 145.0f){
    return 10;
  }else{
    return 11;
  }
}


//returns a bit-pattern specifying which TIP and TIGRESS events
//passed the timing gates specified in common.h
//bits 0 to (MAXNUMTIPHIT-1): TIP hit indices
//bits MAXNUMTIPHIT to MAXNUMTIPHIT+MAXNUMTIGHIT: TIGRESS hit indices
//minTigHit: minimum number of TIGRESS hits in the event (<=1 means no TIG-TIG gate)
//minTipHit: minimum number of TIP hits in the event (<=1 means no TIP-TIP gate)
uint64_t passesTimeGate(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit){

  //Defining Pointers
  TTigressHit *add_hit, *add_hit2;

  bool goodTipTipTime = false;
  bool goodTigTigTime = false;
  bool goodTipTigTime = false;
  uint64_t pass = 0;

  //evaluate timing conditions
  if(tigress){

    if(!tip && (minTipHit>0)){
      return 0;
    }
    if(tip && (tip->GetMultiplicity()>MAXNUMTIPHIT || tip->GetMultiplicity()<minTipHit)){
      return 0;
    }
    if(tigress->GetAddbackMultiplicity()>MAXNUMTIGHIT || tigress->GetAddbackMultiplicity()<minTigHit){
      return 0;
    }
    if((minTigHit==0) && (minTipHit==0)){
      //no timing condition, pass all valid TIP and TIGRESS hits
      if(tip){
        for(int i=0;i<tip->GetMultiplicity();i++){
          if(tip->GetTipHit(i)->GetKValue() != noPileupKValue){
            //pileup
            continue;
          }
          if(tip->GetTipHit(i)->GetPID()<=0){
            //failed fit
            continue;
          }
          pass |= (1ULL<<i);
          goodTipTipTime = true;
        }
      }
      for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
        add_hit = tigress->GetAddbackHit(i);
        if(add_hit->GetKValue() != noPileupKValue){
          //pileup
          continue;
        }
        if(add_hit->BGOFired() || add_hit->GetEnergy() <= 15){
          continue;
        }
        pass |= (1ULL<<(MAXNUMTIPHIT+i));
        goodTigTigTime = true;
      }
      if(goodTipTipTime && goodTigTigTime){
        goodTipTigTime = true;
      }
    }else{
      //evaluate timing conditions

      //evaluate TIP timing
      if(tip){
        if(tip->GetMultiplicity()>=minTipHit){
          if(minTipHit>= 2 && tip->GetMultiplicity()>=2){
            
            //TIP-TIP timing
            for(int i=0;i<tip->GetMultiplicity();i++){
              if(tip->GetTipHit(i)->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }
              if(tip->GetTipHit(i)->GetPID()<=0){
                //failed fit
                continue;
              }
              for(int j=i+1;j<tip->GetMultiplicity();j++){
                if(tip->GetTipHit(j)->GetKValue() != noPileupKValue){
                  //pileup
                  continue;
                }
                if(tip->GetTipHit(j)->GetPID()<=0){
                  //failed fit
                  continue;
                }
                double_t fitT1 = getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger);
                double_t fitT2 = getTipFitTime(tip->GetTipHit(j),tip_waveform_pretrigger);
                if((fitT1 != 0.)&&(fitT2 != 0.)){
                  double_t fitTDiff = fitT1 - fitT2;
                  if(gate1D(fitTDiff,tiptipTGate[0],tiptipTGate[1])){
                    goodTipTipTime = true;
                    pass |= (1ULL<<i);
                    pass |= (1ULL<<j);
                  }
                }
              }
            }
          }else{
            //pass all valid TIP hits
            for(int i=0;i<tip->GetMultiplicity();i++){
              if(tip->GetTipHit(i)->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }
              if(tip->GetTipHit(i)->GetPID()<=0){
                //failed fit
                continue;
              }
              pass |= (1ULL<<i);
              goodTipTipTime = true;
            }
          }
        }
      }

      //evaluate TIGRESS timing
      if(goodTipTipTime || (minTipHit == 0)){
        if(tigress->GetAddbackMultiplicity()>=minTigHit){
          if(minTigHit>=2 && tigress->GetAddbackMultiplicity()>=2){
            for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
              //TIGRESS-TIGRESS addback timing
              add_hit = tigress->GetAddbackHit(i);
              /*if(add_hit->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }*/
              //cout << "energy: " << add_hit->GetEnergy() << ", array num: " << add_hit->GetArrayNumber() << ", address: " << add_hit->GetAddress() << endl;
              if(!add_hit->BGOFired() && add_hit->GetEnergy() > 15){
                for(int j=i+1;j<tigress->GetAddbackMultiplicity();j++){
                  add_hit2 = tigress->GetAddbackHit(j);
                  /*if(add_hit2->GetKValue() != noPileupKValue){
                    //pileup
                    continue;
                  }*/
                  if(!add_hit2->BGOFired() && add_hit2->GetEnergy() > 15){
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
          }else{
            //pass all valid TIGRESS hits
            for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
              add_hit = tigress->GetAddbackHit(i);
              if(add_hit->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }
              if(add_hit->BGOFired() || add_hit->GetEnergy() <= 15){
                continue;
              }
              pass |= (1ULL<<(MAXNUMTIPHIT+i));
              goodTigTigTime = true;
            }
          }
        }
      }

      //evaluate TIP-TIGRESS timing
      if(goodTigTigTime && (minTipHit > 0) && tip){

        bool tigHitTipCoinc[MAXNUMTIGHIT];
        for(int i=0;i<MAXNUMTIGHIT;i++){
          tigHitTipCoinc[i] = false;
        }
        
        for(int i=0;i<tip->GetMultiplicity();i++){
          if(pass&(1ULL<<i)){
            bool tipHitTigCoinc = false;
            for(int j=0;j<tigress->GetAddbackMultiplicity();j++){
              if(pass&(1ULL<<(MAXNUMTIPHIT+j))){
                //K, BGO and energy conditions previously evaluated
                double_t tiptigTDiff = getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger) - tigress->GetAddbackHit(j)->GetTime();
                if(gate1D(tiptigTDiff,tiptigTGate[0],tiptigTGate[1])){
                  tipHitTigCoinc = true;
                  tigHitTipCoinc[j] = true;
                }
              }
            }
            if(tipHitTigCoinc == false){
              pass &= ~(1ULL<<i); //unset
            }
          }
        }

        for(int i=0;i<MAXNUMTIGHIT;i++){
          if(tigHitTipCoinc[i] == false){
            pass &= ~(1ULL<<(MAXNUMTIPHIT+i)); //unset
          }
        }

        //count the number of passed hits
        uint8_t numPassedTip = 0;
        uint8_t numPassedTig = 0;
        for(int i=0;i<MAXNUMTIPHIT;i++){
          if(pass & (1ULL<<i)){
            numPassedTip++;
          }
        }
        for(int i=0;i<MAXNUMTIGHIT;i++){
          if(pass & (1ULL<<(MAXNUMTIPHIT+i))){
            numPassedTig++;
          }
        }
        if((numPassedTip>=minTipHit)&&(numPassedTig>=minTigHit)){
          goodTipTigTime = true;
        }
      }

      if(goodTipTigTime){
        pass |= (1ULL<<TIPTIGFLAG);
      }
      if(goodTigTigTime){
        pass |= (1ULL<<TIGTIGFLAG);
      }
      if(goodTipTipTime){
        pass |= (1ULL<<TIPTIPFLAG);
      }
    }

  }
  
  /*if(goodTipTipTime)
  if(goodTigTigTime)
  if(!goodTipTigTime){
    //report failed events
    cout << "Event failed timing:" << endl;
    for(int i=0;i<tip->GetMultiplicity();i++){
      cout << " TIP hit " << i << ", t:" << tip->GetTipHit(i)->GetTime() << ", K: " << tip->GetTipHit(i)->GetKValue() << ", PID: " << tip->GetTipHit(i)->GetPID() << endl;
    }
    for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetAddbackMultiplicity(); tigHitIndAB++){
      cout << " TIGRESS AB hit " << tigHitIndAB << ", t: " << tigress->GetAddbackHit(tigHitIndAB)->GetTime() << ", K: " 
      << tigress->GetAddbackHit(tigHitIndAB)->GetKValue() << ", BGO: " << tigress->GetAddbackHit(tigHitIndAB)->BGOFired() << ", E: " << tigress->GetAddbackHit(tigHitIndAB)->GetEnergy() << endl;
    }
    if(goodTipTipTime)
      cout << "TIP-TIP timing passed" << endl;
    if(goodTigTigTime)
      cout << "TIG-TIG timing passed" << endl;
  }*/
  

  return pass;
  
}
