//File contains common functions which all other sort codes
//have access to

#define common_cxx
#include "../include/common.h"
#include "position_arrays.cxx"

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
        protonRingCut[i]->SetPoint(1,2.9,125);
        protonRingCut[i]->SetPoint(2,14,143);
        protonRingCut[i]->SetPoint(3,14,82.9);
        protonRingCut[i]->SetPoint(4,7.7,43.5);
        protonRingCut[i]->SetPoint(5,5.6,29.3);
        protonRingCut[i]->SetPoint(6,3.7,27.7);
        protonRingCut[i]->SetPoint(7,2,51);

        alphaRingCut[i]->SetPoint(0,2.2,24);
        alphaRingCut[i]->SetPoint(1,5,24);
        alphaRingCut[i]->SetPoint(2,20,37);
        alphaRingCut[i]->SetPoint(3,20,13);
        alphaRingCut[i]->SetPoint(4,8,2);
        alphaRingCut[i]->SetPoint(5,3.7,2);
        alphaRingCut[i]->SetPoint(6,2,11);
        alphaRingCut[i]->SetPoint(7,2.2,24);
        break;
      case 8:
        //ring 8
        protonRingCut[i]->SetPoint(0,1,37);
        protonRingCut[i]->SetPoint(1,3.2,110);
        protonRingCut[i]->SetPoint(2,14,130);
        protonRingCut[i]->SetPoint(3,14,83.6);
        protonRingCut[i]->SetPoint(4,7.5,40);
        protonRingCut[i]->SetPoint(5,7,29.5);
        protonRingCut[i]->SetPoint(6,3.6,24.1);
        protonRingCut[i]->SetPoint(7,1,37);

        alphaRingCut[i]->SetPoint(0,1.2,22.8);
        alphaRingCut[i]->SetPoint(1,4.9,23.7);
        alphaRingCut[i]->SetPoint(2,19,27);
        alphaRingCut[i]->SetPoint(3,19,20);
        alphaRingCut[i]->SetPoint(4,9,3);
        alphaRingCut[i]->SetPoint(5,4.2,3);
        alphaRingCut[i]->SetPoint(6,2.7,6.5);
        alphaRingCut[i]->SetPoint(7,1.2,22.8);
        break;
      case 7:
        //ring 7
        protonRingCut[i]->SetPoint(0,2.5,33.5);
        protonRingCut[i]->SetPoint(1,2.7,121.5);
        protonRingCut[i]->SetPoint(2,20,144);
        protonRingCut[i]->SetPoint(3,20,76);
        protonRingCut[i]->SetPoint(4,9.8,55);
        protonRingCut[i]->SetPoint(5,7.6,45);
        protonRingCut[i]->SetPoint(6,6.75,37);
        protonRingCut[i]->SetPoint(7,2.5,33.5);

        alphaRingCut[i]->SetPoint(0,1.9,26);
        alphaRingCut[i]->SetPoint(1,26,46);
        alphaRingCut[i]->SetPoint(2,27,15);
        alphaRingCut[i]->SetPoint(3,8,2);
        alphaRingCut[i]->SetPoint(4,3,2);
        alphaRingCut[i]->SetPoint(5,1.5,6);
        alphaRingCut[i]->SetPoint(6,1.2,13);
        alphaRingCut[i]->SetPoint(7,1.9,26);
        break;
      case 6:
        //ring 6
        protonRingCut[i]->SetPoint(0,1.8,31);
        protonRingCut[i]->SetPoint(1,2.5,105);
        protonRingCut[i]->SetPoint(2,24,126);
        protonRingCut[i]->SetPoint(3,24,78);
        protonRingCut[i]->SetPoint(4,12,64);
        protonRingCut[i]->SetPoint(5,9,46);
        protonRingCut[i]->SetPoint(6,7,31);
        protonRingCut[i]->SetPoint(7,1.8,31);

        alphaRingCut[i]->SetPoint(0,1.8,22);
        alphaRingCut[i]->SetPoint(1,7.5,29);
        alphaRingCut[i]->SetPoint(2,32,32);
        alphaRingCut[i]->SetPoint(3,33,11);
        alphaRingCut[i]->SetPoint(4,10,2);
        alphaRingCut[i]->SetPoint(5,3,2);
        alphaRingCut[i]->SetPoint(6,1.5,7);
        alphaRingCut[i]->SetPoint(7,1.8,22);
        break;
      case 5:
        //ring 5
        protonRingCut[i]->SetPoint(0,1.2,36);
        protonRingCut[i]->SetPoint(1,2.9,112);
        protonRingCut[i]->SetPoint(2,20,124);
        protonRingCut[i]->SetPoint(3,20,84);
        protonRingCut[i]->SetPoint(4,11.4,60);
        protonRingCut[i]->SetPoint(5,9,49.5);
        protonRingCut[i]->SetPoint(6,8.3,44.8);
        protonRingCut[i]->SetPoint(7,1.2,36);

        alphaRingCut[i]->SetPoint(0,1.9,30);
        alphaRingCut[i]->SetPoint(1,5.8,38);
        alphaRingCut[i]->SetPoint(2,21.5,61);
        alphaRingCut[i]->SetPoint(3,21,12);
        alphaRingCut[i]->SetPoint(4,9.5,2);
        alphaRingCut[i]->SetPoint(5,2.6,2);
        alphaRingCut[i]->SetPoint(6,1,11);
        alphaRingCut[i]->SetPoint(7,1.9,30);
        break;
      case 4:
        //ring 4
        protonRingCut[i]->SetPoint(0,1.8,31.5);
        protonRingCut[i]->SetPoint(1,3,120);
        protonRingCut[i]->SetPoint(2,17.5,135);
        protonRingCut[i]->SetPoint(3,17.5,83.5);
        protonRingCut[i]->SetPoint(4,11.8,57);
        protonRingCut[i]->SetPoint(5,9,44.8);
        protonRingCut[i]->SetPoint(6,6.2,40.5);
        protonRingCut[i]->SetPoint(7,1.8,31.5);

        alphaRingCut[i]->SetPoint(0,2.6,28.2);
        alphaRingCut[i]->SetPoint(1,8,41);
        alphaRingCut[i]->SetPoint(2,23,60);
        alphaRingCut[i]->SetPoint(3,23,25);
        alphaRingCut[i]->SetPoint(4,10.6,9.5);
        alphaRingCut[i]->SetPoint(5,8.5,2);
        alphaRingCut[i]->SetPoint(6,2.3,2);
        alphaRingCut[i]->SetPoint(7,2.6,28.2);
        break;
      case 3:
        //ring 3
        protonRingCut[i]->SetPoint(0,1.5,44);
        protonRingCut[i]->SetPoint(1,3.2,130);
        protonRingCut[i]->SetPoint(2,12,146);
        protonRingCut[i]->SetPoint(3,20,150);
        protonRingCut[i]->SetPoint(4,21.2,87);
        protonRingCut[i]->SetPoint(5,12.8,63.5);
        protonRingCut[i]->SetPoint(6,11.6,50);
        protonRingCut[i]->SetPoint(7,1.5,45);

        alphaRingCut[i]->SetPoint(0,2.9,42);
        alphaRingCut[i]->SetPoint(1,10,47);
        alphaRingCut[i]->SetPoint(2,28,59);
        alphaRingCut[i]->SetPoint(3,28,26);
        alphaRingCut[i]->SetPoint(4,12.8,11.5);
        alphaRingCut[i]->SetPoint(5,9.6,2);
        alphaRingCut[i]->SetPoint(6,1.5,2);
        alphaRingCut[i]->SetPoint(7,2.9,42);
        break;
      case 2:
        //ring 2
        protonRingCut[i]->SetPoint(0,1.5,43);
        protonRingCut[i]->SetPoint(1,3,134);
        protonRingCut[i]->SetPoint(2,23,146);
        protonRingCut[i]->SetPoint(3,25,98);
        protonRingCut[i]->SetPoint(4,16,71);
        protonRingCut[i]->SetPoint(5,13.5,53);
        protonRingCut[i]->SetPoint(6,8,45);
        protonRingCut[i]->SetPoint(7,1.5,43);

        alphaRingCut[i]->SetPoint(0,2.2,41);
        alphaRingCut[i]->SetPoint(1,11,44);
        alphaRingCut[i]->SetPoint(2,30,53);
        alphaRingCut[i]->SetPoint(3,31,31);
        alphaRingCut[i]->SetPoint(4,14.5,2);
        alphaRingCut[i]->SetPoint(5,4,2);
        alphaRingCut[i]->SetPoint(6,1,4);
        alphaRingCut[i]->SetPoint(7,2.2,41);
        break;
      case 1:
        //ring 1
        protonRingCut[i]->SetPoint(0,3,103);
        protonRingCut[i]->SetPoint(1,8,126);
        protonRingCut[i]->SetPoint(2,30,127);
        protonRingCut[i]->SetPoint(3,35,91);
        protonRingCut[i]->SetPoint(4,22,55);
        protonRingCut[i]->SetPoint(5,12,49);
        protonRingCut[i]->SetPoint(6,2.1,46.5);
        protonRingCut[i]->SetPoint(7,3,103);

        alphaRingCut[i]->SetPoint(0,2.7,44);
        alphaRingCut[i]->SetPoint(1,22,54);
        alphaRingCut[i]->SetPoint(2,43,65);
        alphaRingCut[i]->SetPoint(3,43,29);
        alphaRingCut[i]->SetPoint(4,21,13);
        alphaRingCut[i]->SetPoint(5,15.6,2);
        alphaRingCut[i]->SetPoint(6,2.4,2);
        alphaRingCut[i]->SetPoint(7,2.7,44);
        break;
      case 0:
      default:
        //ring 0
        protonRingCut[i]->SetPoint(0,3,50);
        protonRingCut[i]->SetPoint(1,3.3,127);
        protonRingCut[i]->SetPoint(2,17,158);
        protonRingCut[i]->SetPoint(3,26,138);
        protonRingCut[i]->SetPoint(4,26,73);
        protonRingCut[i]->SetPoint(5,12,57);
        protonRingCut[i]->SetPoint(6,7.8,50);
        protonRingCut[i]->SetPoint(7,3,50);

        alphaRingCut[i]->SetPoint(0,3,40);
        alphaRingCut[i]->SetPoint(1,8,48);
        alphaRingCut[i]->SetPoint(2,18,62);
        alphaRingCut[i]->SetPoint(3,32,76);
        alphaRingCut[i]->SetPoint(4,30,23);
        alphaRingCut[i]->SetPoint(5,16.1,4);
        alphaRingCut[i]->SetPoint(6,2.6,3);
        alphaRingCut[i]->SetPoint(7,3,40);
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

TVector3 getTigVector(uint8_t core, uint8_t seg){
  TVector3 hitPos(0,0,0);
  if(core > 63){
    printf("WARNING: bad core value (%u)\n",core);
  }else if(seg > 8){
    printf("WARNING: bad segment value (%u)\n",seg);
  }else{
    switch(core % 4){
      case 3:
        //white core
        hitPos.SetXYZ(GeWhitePositionBack[(core/4) + 1][seg][0],GeWhitePositionBack[(core/4) + 1][seg][1],GeWhitePositionBack[(core/4) + 1][seg][2]);
        break;
      case 2:
        //red core
        hitPos.SetXYZ(GeRedPositionBack[(core/4) + 1][seg][0],GeRedPositionBack[(core/4) + 1][seg][1],GeRedPositionBack[(core/4) + 1][seg][2]);
        break;
      case 1:
        //green core
        hitPos.SetXYZ(GeGreenPositionBack[(core/4) + 1][seg][0],GeGreenPositionBack[(core/4) + 1][seg][1],GeGreenPositionBack[(core/4) + 1][seg][2]);
        break;
      case 0:
      default:
        //blue core
        hitPos.SetXYZ(GeBluePositionBack[(core/4) + 1][seg][0],GeBluePositionBack[(core/4) + 1][seg][1],GeBluePositionBack[(core/4) + 1][seg][2]);
        break;
    }
  }
  return hitPos;
}

double_t TigGetDoppler(double beta, double eTig, uint8_t core, uint8_t seg, TVector3* vec = nullptr){
  if(vec == nullptr) {
    cout << "ERROR: need to define momentum vector for Doppler correction!" << endl;
    exit(-1);
  }
  if(seg > 8){
    cout << "ERROR: invalid TIGRESS segment (" << seg << ")." << endl;
    exit(-1);
  }
  double tmp   = 0;
  double gamma = 1.0 / (sqrt(1.0 - pow(beta, 2)));
  TVector3 hitPos = getTigVector(core,seg);
  
  tmp = eTig * gamma * (1 - beta * TMath::Cos(hitPos.Angle(*vec)));
  return tmp;
}

//for when we want to manually specify hit properties (eg. when summing non-addback data)
double_t getEDoppFusEvapManualBeta(double eTig, uint8_t core, uint8_t seg, uint8_t numCsIHits, Double_t beta, csi_hit *tip_hits, PIDGates *gates){
  double_t resM = compoundM_AMU*AMU; //residual mass prior to particle evaporation
  TVector3 p_compound(0,0,resM*beta/(1.0-beta*beta));
  //cout << "p_compound: " << p_compound.X() << " " << p_compound.Y() << " " << p_compound.Z() << endl;
  TVector3 p_part;
  double_t resBeta = 0.;
  if(numCsIHits<=MAXNUMTIPHIT){
    for(int tipHitInd=0;tipHitInd<numCsIHits;tipHitInd++){
      switch(getParticleTypePID(tip_hits[tipHitInd].PID,tip_hits[tipHitInd].energy,tip_hits[tipHitInd].detNum,gates)){ //see common.cxx
        case 4:
          p_part.SetXYZ(TipCsIPosition[tip_hits[tipHitInd].detNum-1][0],TipCsIPosition[tip_hits[tipHitInd].detNum-1][1],TipCsIPosition[tip_hits[tipHitInd].detNum-1][2]);
          //cout << "p_alpha " << tipHitInd << ": " << p_part.X() << " " << p_part.Y() << " " << p_part.Z() << endl;
          p_part.SetMag(sqrt(2.0*4.0015*AMU*tip_hits[tipHitInd].energy));
          //p_part.SetMag(sqrt(2.0*4.0015*AMU*10)); //TEMPORARY - ASSUME 10 MeV
          p_compound -= p_part;
          resM -= 2.0*4.0015*AMU; //technically not true due to binding energy!
          break;
        case 1:
          p_part.SetXYZ(TipCsIPosition[tip_hits[tipHitInd].detNum-1][0],TipCsIPosition[tip_hits[tipHitInd].detNum-1][1],TipCsIPosition[tip_hits[tipHitInd].detNum-1][2]);
          //cout << "p_proton " << tipHitInd << ": " << p_part.X() << " " << p_part.Y() << " " << p_part.Z() << endl;
          p_part.SetMag(sqrt(2.0*938.272*tip_hits[tipHitInd].energy));
          //p_part.SetMag(sqrt(2.0*938.272*10)); //TEMPORARY - ASSUME 10 MeV
          p_compound -= p_part;
          resM -= 2.0*938.272; //technically not true due to binding energy!
          break;
        case 0:
        default:
          break;
      }
    }
  }
  //cout << "p_res: " << p_compound.X() << " " << p_compound.Y() << " " << p_compound.Z() << endl;
  resBeta = sqrt(1.0 - 1.0/(1.0 + pow(p_compound.Mag()/resM,2)));
  //cout << "resBeta: " << resBeta << ", E: " << TigGetDoppler(resBeta,add_hit->energy,add_hit->core,add_hit->seg,&p_compound) << endl;
  return TigGetDoppler(resBeta,eTig,core,seg,&p_compound);
}

//for when we want to manually specify hit properties (eg. when summing non-addback data)
double_t getEDoppFusEvapManual(double eTig, uint8_t core, uint8_t seg, uint8_t numCsIHits, csi_hit *tip_hits, PIDGates *gates){
  return getEDoppFusEvapManualBeta(eTig,core,seg,numCsIHits,betaCompound,tip_hits,gates);
}

//implementation of getEDoppFusEvap for the SMOL data format
double_t getEDoppFusEvapDirect(tig_hit *add_hit, uint8_t numCsIHits, csi_hit *tip_hits, PIDGates *gates){
  return getEDoppFusEvapManual(add_hit->energy,add_hit->core,add_hit->seg,numCsIHits,tip_hits,gates);
}

//implementation of getEDoppFusEvap for the SMOL data format
double_t getEDoppFusEvapDirectBeta(tig_hit *add_hit, uint8_t numCsIHits, Double_t beta, csi_hit *tip_hits, PIDGates *gates){
  return getEDoppFusEvapManualBeta(add_hit->energy,add_hit->core,add_hit->seg,numCsIHits,beta,tip_hits,gates);
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


//Get the segment 'ring'
//there are actually 8 unique segment angles per core, clustered
//over 'rings' each with 4 angles covering a ~2.5 deg range
//there are also core-only hits without any assigned segment
Int_t getTIGRESSSegmentRing(const float theta){
  if((theta >= 30.0f)&&(theta <= 37.0f)){
    return 0; //34.325 deg (low 32.95, high 35.55)
  }else if((theta >= 40.0f)&&(theta <= 43.0f)){
    return 1; //41.283 deg (low 40.55, high 42.15, missing 1 angle)
  }else if((theta >= 48.0f)&&(theta <= 52.0f)){
    return 2; //49.875 deg (low 49.15, high 50.55)
  }else if((theta >= 55.0f)&&(theta <= 59.0f)){
    return 3; //57.075 deg (low 56.05, high 58.25)
  }else if((theta >= 77.0f)&&(theta <= 80.0f)){
    return 4; //78.50 deg (low 77.75, high 79.35)
  }else if((theta >= 84.0f)&&(theta <= 88.0f)){
    return 5; //85.683 deg (low 85.35, high 86.05, missing 1 angle)
  }else if((theta >= 92.0f)&&(theta <= 96.0f)){
    return 6; //94.317 deg (low 93.95, high 94.65, missing 1 angle)
  }else if((theta >= 100.0f)&&(theta <= 103.0f)){
    return 7; //101.50 deg (low 100.65, high 102.25)
  }else if((theta >= 121.0f)&&(theta <= 125.0f)){
    return 8; //122.920 deg (low 121.75, high 123.95)
  }else if((theta >= 128.0f)&&(theta <= 132.0f)){
    return 9; //130.120 deg (low 129.45, high 130.85)
  }else if((theta >= 137.0f)&&(theta <= 141.0f)){
    return 10; //138.72 deg (low 137.85, high 139.45, missing 1 angle)
  }else if((theta >= 143.0f)&&(theta <= 148.0f)){
    return 11; //145.680 deg (low 144.45, high 147.05)
  }else{
    //cout << "Invalid TIGRESS angle: " << theta << endl;
    return 12; //no segment
  }
}

Int_t getTIGRESSPhiRing(const float theta, const float phi){
  Int_t ring = getTIGRESSSegmentRing(theta);
  if(ring < 12){
    if((phi > -180)&&(phi < -167)){
      return 0; 
    }else if((phi > -165)&&(phi < -160)){
      return 1; 
    }else if((phi > -155)&&(phi < -148)){
      return 2; 
    }else if((phi > -148)&&(phi < -134)){
      return 3;
    }else if((phi > -127)&&(phi < -122)){
      return 4;
    }else if((phi > -118)&&(phi < -114)){
      return 5;
    }else if((phi > -110)&&(phi < -106)){
      return 6;
    }else if((phi > -103)&&(phi < -100)){
      return 7;
    }else if((phi > -91)&&(phi < -77)){
      return 8;
    }else if((phi > -77)&&(phi < -70)){
      return 9;
    }else if((phi > -65)&&(phi < -58)){
      return 10;
    }else if((phi > -58)&&(phi < -45)){
      return 11;
    }else if((phi > -36)&&(phi < -32)){
      return 12;
    }else if((phi > -28)&&(phi < -25)){
      return 13;
    }else if((phi > -20)&&(phi < -16)){
      return 14;
    }else if((phi > -13)&&(phi < -8)){
      return 15;
    }else if((phi > 0)&&(phi < 13)){
      return 16;
    }else if((phi > 13)&&(phi < 22)){
      return 17;
    }else if((phi > 22)&&(phi < 32)){
      return 18;
    }else if((phi > 32)&&(phi < 46)){
      return 19;
    }else if((phi > 53)&&(phi < 59)){
      return 20;
    }else if((phi > 61)&&(phi < 66)){
      return 21;
    }else if((phi > 69)&&(phi < 74)){
      return 22;
    }else if((phi > 77)&&(phi < 82)){
      return 23;
    }else if((phi > 89)&&(phi < 103)){
      return 24;
    }else if((phi > 103)&&(phi < 114)){
      return 25;
    }else if((phi > 114)&&(phi < 122)){
      return 26;
    }else if((phi > 122)&&(phi < 136)){
      return 27;
    }else if((phi > 143)&&(phi < 149)){ //estimated, no detector present here during S2182
      return 28;
    }else if((phi > 151)&&(phi < 156)){ //estimated, no detector present here during S2182
      return 29;
    }else if((phi > 159)&&(phi < 164)){ //estimated, no detector present here during S2182
      return 30;
    }else if((phi > 167)&&(phi < 172)){ //estimated, no detector present here during S2182
      return 31;
    }else{
      //cout << "Invalid TIGRESS angle: " << theta << endl;
      return 32; //no phi ring
    }
  }else{
    //cout << "Invalid TIGRESS angle: " << theta << endl;
    return 32; //no phi ring
  }
}


uint64_t passesTimeGateNoAB(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit){

  //Defining Pointers
  TTigressHit *tig_hit, *tig_hit2;

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
    if(tigress->GetMultiplicity()>MAXNUMTIGHIT || tigress->GetMultiplicity()<minTigHit){
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
      for(int i=0;i<tigress->GetMultiplicity();i++){
        tig_hit = tigress->GetTigressHit(i);
        if(tig_hit->GetKValue() != noPileupKValue){
          //pileup
          continue;
        }
        if(tig_hit->BGOFired() || tig_hit->GetEnergy() <= 15){
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
      int mult = tigress->GetMultiplicity();
      if(goodTipTipTime || (minTipHit == 0)){
        if(mult>=minTigHit){
          if(minTigHit>=2 && mult>=2){
            for(int i=0;i<tigress->GetMultiplicity();i++){
              //TIGRESS-TIGRESS non-addback timing
              tig_hit = tigress->GetTigressHit(i);
              /*if(tig_hit->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }*/
              //cout << "energy: " << tig_hit->GetEnergy() << ", array num: " << tig_hit->GetArrayNumber() << ", address: " << tig_hit->GetAddress() << endl;
              if(!tig_hit->BGOFired() && tig_hit->GetEnergy() > 15){
                for(int j=i+1;j<tigress->GetMultiplicity();j++){
                  tig_hit2 = tigress->GetTigressHit(j);
                  /*if(tig_hit2->GetKValue() != noPileupKValue){
                    //pileup
                    continue;
                  }*/
                  if(!tig_hit2->BGOFired() && tig_hit2->GetEnergy() > 15){
                    double tigtigTDiff = tig_hit->GetTime() - tig_hit2->GetTime();
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
            for(int i=0;i<tigress->GetMultiplicity();i++){
              tig_hit = tigress->GetTigressHit(i);
              if(tig_hit->GetKValue() != noPileupKValue){
                //pileup
                continue;
              }
              if(tig_hit->BGOFired() || tig_hit->GetEnergy() <= 15){
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
            for(int j=0;j<tigress->GetMultiplicity();j++){
              if(pass&(1ULL<<(MAXNUMTIPHIT+j))){
                //K, BGO and energy conditions previously evaluated
                double_t tiptigTDiff = getTipFitTime(tip->GetTipHit(i),tip_waveform_pretrigger) - tigress->GetTigressHit(j)->GetTime();
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
    for(int tigHitIndAB = 0; tigHitIndAB < tigress->GetMultiplicity(); tigHitIndAB++){
      cout << " TIGRESS hit " << tigHitIndAB << ", t: " << tigress->GetTigressHit(tigHitIndAB)->GetTime() << ", K: " 
      << tigress->GetTigressHit(tigHitIndAB)->GetKValue() << ", BGO: " << tigress->GetTigressHit(tigHitIndAB)->BGOFired() << ", E: " << tigress->GetTigressHit(tigHitIndAB)->GetEnergy() << endl;
    }
    if(goodTipTipTime)
      cout << "TIP-TIP timing passed" << endl;
    if(goodTigTigTime)
      cout << "TIG-TIG timing passed" << endl;
  }*/
  

  return pass;
  
}


//returns a bit-pattern specifying which TIP and TIGRESS events
//passed the timing gates specified in common.h
//bits 0 to (MAXNUMTIPHIT-1): TIP hit indices
//bits MAXNUMTIPHIT to MAXNUMTIPHIT+MAXNUMTIGHIT: TIGRESS hit indices
//minTigHit: minimum number of TIGRESS hits in the event (<=1 means no TIG-TIG gate)
//minTipHit: minimum number of TIP hits in the event (<=1 means no TIP-TIP gate)
//noAddback: whether to use addback or non-addback hits
uint64_t passesTimeGateAB(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit){

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
      int mult = tigress->GetAddbackMultiplicity();
      if(goodTipTipTime || (minTipHit == 0)){
        if(mult>=minTigHit){
          if(minTigHit>=2 && mult>=2){
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

uint64_t passesTimeGate(TTigress *tigress, TTip *tip, const uint8_t minTigHit, const uint8_t minTipHit){
  return passesTimeGateAB(tigress,tip,minTigHit,minTipHit);
}
