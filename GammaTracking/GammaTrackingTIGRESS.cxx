#include "GammaTrackingTIGRESS.h" //define all global variables here!

TRandom3 *randGen;

Int_t getNumAngleBins(Int_t rInd, Double_t rScaleFac, Double_t scaleFac){ 
  return 1 + (Int_t)(pow((rInd/((VOXEL_BINS_R*rScaleFac) - 1.0)),2.0)*((VOXEL_BINS_ANGLE_MAX*scaleFac)-1)); //the number of angle bins in the map depends on r
} 

Int_t getNumAngleBins(Int_t rInd, Double_t scaleFac){ 
  return getNumAngleBins(rInd,scaleFac,scaleFac); 
}

//function for calculation of ordering parameters
//this enforces the single segment hit condition assumed by the mapping process, using SEGMENT_ENERGY_THRESHOLD
double calc_ordering(TTigressHit * tigress_hit, const Int_t i, const Int_t parameterNum) {

  //cout << "array position: " << tigress_hit->GetArrayNumber() << endl;

  //lists of adjacent segments in the TIGRESS array (zero-indexed)
  const Int_t phiAdjSeg1[8] = {3,0,1,2,7,4,5,6};
  const Int_t phiAdjSeg2[8] = {1,2,3,0,5,6,7,4};
  const Int_t    zAdjSeg[8] = {4,5,6,7,0,1,2,3};
  
  const TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(i);
  const Int_t numSamples = segment_hit.GetWaveform()->size();

  if(segment_hit.GetCharge() < SEGMENT_ENERGY_THRESHOLD){
    return BAD_RETURN;
  }
  if(segment_hit.GetCharge() > MAX_ENERGY_SINGLE_INTERACTION){
    return BAD_RETURN; //energy too high - event likely to result from multiple hits
  }

  const Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed
  //cout << "segment " << segNum << " energy: " << segment_hit.GetEnergy() << endl;

  const std::vector<Short_t> *segwf;
  segwf = segment_hit.GetWaveform();
  if(segwf->size() != numSamples){
    //cout << "Mismatched waveform sizes." << endl;
    return BAD_RETURN;
  }

  if(parameterNum==0){
    //construct rho, the ordering parameter for the radius
    //see Eq. 4 of NIM A 729 (2013) 198-206
    
    const std::vector<Short_t> *segwf2, *segwf3, *segwf4;
    bool found1 = false;
    bool found2 = false;
    bool found3 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if(j!=i){
        if(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD){
          return BAD_RETURN;
        }
        //cout << "energy: " << tigress_hit->GetSegmentHit(j).GetEnergy() << ", charge: " << tigress_hit->GetSegmentHit(j).GetCharge() << endl;
        if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
          if(found1==false){
            found1=true;
            segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }else{
            return BAD_RETURN;
          }
        }
        if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
          if(found2==false){
            found2=true;
            segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }else{
            return BAD_RETURN;
          }
        }
        if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
          if(found3==false){
            found3=true;
            segwf4 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }else{
            return BAD_RETURN;
          }
        }
      }
    }
    if((!found1)||(!found2)||(!found3)){
      //cout << "Cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
      return BAD_RETURN;
    }else if((segwf2->size() != numSamples)||(segwf3->size() != numSamples)||(segwf4->size() != numSamples)){
      //cout << "Mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    double dno = 0.; //placeholder for denominator
    double sampleAvg = 0.;
    double term2 = 0.;
    for(int j=1;j<numSamples-1;j++){
      sampleAvg += (j)*(segwf->at(j+1) - segwf->at(j-1))/2.0;
      term2 += segwf2->at(j) + segwf3->at(j);
      //term2 += segwf2->at(j) + segwf3->at(j) + segwf4->at(j);
      dno += (segwf->at(j+1) - segwf->at(j-1))/2.0;
    }
    sampleAvg /= dno;
    //cout << "sampleAvg: " << sampleAvg << endl;
    term2 /= dno*2.0;
    double rho = 0.;
    for(int j=1;j<numSamples-1;j++){
      rho += pow(j - sampleAvg,3.0)*(segwf->at(j+1) - segwf->at(j-1))/2.0;
    }
    rho /= dno;
    //rho = rho - term2; //may want to adjust weight of term2, Li uses x400
    //rho = -1.*term2;
    if((dno==0.)||(rho!=rho)){
    //if(rho!=rho){
      //cout << "Cannot compute rho parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    //cout << "term2: " << term2 << ", rho: " << rho << endl;
    return rho;
  }else if(parameterNum==1){
    //contruct phi, the ordering parameter for the angle
    //see Eq. 3 of NIM A 729 (2013) 198-206
    const std::vector<Short_t> *segwf2, *segwf3;
    bool found1 = false;
    bool found2 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if(j!=i){
        if(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD){
          return BAD_RETURN;
        }
        if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
          if(found1==false){
            found1=true;
            segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }else{
            return BAD_RETURN;
          }
        }
        if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
          if(found2==false){
            found2=true;
            segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
          }else{
            return BAD_RETURN;
          }
        }
      }
    }
    if((!found1)||(!found2)){
      //cout << "Cannot get neighbouring segment wavefoms to compute phi parameter." << endl;
      return BAD_RETURN;
    }else if((segwf2->size() != numSamples)||(segwf3->size() != numSamples)){
      //cout << "Mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    double phi = 0.;
    //use the Li 2018 method here, the Desequelles 2013 method fails when sorting one source with a map from a different source
    double maxVall = -1E30;
    double minVall = 1E30;
    double maxValr = -1E30;
    double minValr = 1E30;
    for(int j=0;j<numSamples ;j++){
      if(segwf2->at(j) > maxVall){
        maxVall = segwf2->at(j);
      }
      if(segwf2->at(j) < minVall){
        minVall = segwf2->at(j);
      }
      if(segwf3->at(j) > maxValr){
        maxValr = segwf3->at(j);
      }
      if(segwf3->at(j) < minValr){
        minValr = segwf3->at(j);
      }
    }
    if((minVall == 1E30)||(maxVall == -1E30)||(minValr == 1E30)||(maxValr == -1E30)){
      //cout << "Cannot find maximum or minimum values for phi." << endl;
      //cout << "vals: " << minVall << " " << maxVall << " " << minValr << " " << maxValr << endl;
      return BAD_RETURN;
    }
    phi = ((maxValr - minValr) - (maxVall - minVall))/((maxVall - minVall) + (maxValr - minValr));
    if((phi!=phi)){
      //cout << "Cannot compute phi parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    //cout << "phi: " << phi << endl;
    return phi;
  }else if(parameterNum==2){
    //contruct zeta, the ordering parameter for the z direction
    //see Eq. 2 of NIM A 729 (2013) 198-206 (modified here)

    double zeta = 0.;
    double dno; //placeholder for denominator
    if(segNum>3){
      //back segment
       const std::vector<Short_t> *segwf2;
      bool found1 = false;
      for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
        if(j!=i){
          if(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD){
            return BAD_RETURN;
          }
          if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
            if(found1==false){
              found1=true;
              segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
            }else{
              return BAD_RETURN;
            }
          }
        }
      }
      if(!found1){
        //cout << "Cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
        return BAD_RETURN;
      }else if(segwf2->size() != numSamples){
        //cout << "Mismatched waveform sizes." << endl;
        return BAD_RETURN;
      }
      double maxVal = -1E30;
      double minVal = 1E30;
      dno = 0.;
      for(int j=1;j<numSamples-1;j++){
        if(segwf2->at(j) > maxVal){
          maxVal = segwf2->at(j);
        }
        if(segwf2->at(j) < minVal){
          minVal = segwf2->at(j);
        }
        dno += (segwf->at(j+1) - segwf->at(j-1))/2.0;
      }
      if((minVal == 1E30)||(maxVal == 0)||(maxVal == minVal)){
        //cout << "Cannot find maximum or minimum values for zeta." << endl;
        return BAD_RETURN;
      }
      zeta = (maxVal - minVal)/dno;
      
      //cout << "max: " << maxVal << ", min: " << minVal << ", dno: " << dno << ", zeta: " << zeta << endl;
      zeta *= -2.; //back segment, reverse sign to make zeta increase with z
      zeta += 1.0;
    }else{
      //front segment
      dno = 0.;
      for(int j=1;j<numSamples-1;j++){
        zeta += (j)*(segwf->at(j+1) - segwf->at(j-1))/2.0;
        dno += (segwf->at(j+1) - segwf->at(j-1))/2.0;
      }
      zeta /= dno;
      //get zeta values centered on 0
      zeta *= 4./numSamples;
      zeta -= 2.0;
    }
    if((dno==0.)||(zeta!=zeta)){
      //cout << "Cannot compute zeta parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    //cout << "zeta: " << zeta << endl;
    return zeta;
  }else{
    return BAD_RETURN;
  }

}

//checks if there is a common hit in the hitpatterns,
//and returns the bit index of the first common hit
int32_t commonHitInHP(int32_t hp1, int32_t hp2, int32_t max_search){
  if(max_search>32){
    cout << "WARNING: trying to search past the bounds of a 32-bit integer!" << endl;
    return -1;
  }
  for(int i=0;i<max_search;i++){
    if((hp1&(1<<i))&&(hp2&(1<<i))){
      return i;
    }
  }

  return -1;
}


//import gamma tracking map from a ROOT file
//the map is responsible for mapping ordering parameters to real spatial coordinates
void GT_import_map(TFile *map_file, GT_map *gt_map){

  char hname[64];

  if(!map_file->IsOpen()) {
    cout << "ERROR: GT_import_map - could not open map file!" << endl;
    exit(-1);
  }

  //read in histograms from map file
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      sprintf(hname,"rMapPos%iSeg%i",l,k);
      if((gt_map->rMap[l*NSEG + k] = (TH1*)map_file->Get(hname))==NULL){
        //cout << "No r coordinate map for segment " << k << endl;
      }
      for(int j = 0; j < VOXEL_BINS_R; j++){
        Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
        Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
        sprintf(hname,"angleMapPos%iSeg%ir%ito%i",l,k,rValMin,rValMax);
        if((gt_map->angleMap[l*NSEG*VOXEL_BINS_R + k*VOXEL_BINS_R + j] = (TH1*)map_file->Get(hname))==NULL){
          //cout << "No angle coordinate map for segment " << k << ", radial bin " << j << endl;
        }
        Int_t numAngleBins = getNumAngleBins(j,1.0);
        for(int i = 0; i < numAngleBins; i++){
          Int_t angleValMin = (Int_t)(i*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
          Int_t angleValMax = (Int_t)((i+1)*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
          sprintf(hname,"zMapPos%iSeg%ir%ito%iangle%ito%i",l,k,rValMin,rValMax,angleValMin,angleValMax);
          if((gt_map->zMap[l*NSEG*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = (TH1*)map_file->Get(hname))==NULL){
            gt_map->zMap[l*NSEG*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]=NULL;
            //cout << "No z coordinate map for segment " << k << ", radial bin " << j << ", angle bin " << i << ", name: " << hname << endl;
          }
        }
      } 
    }
  }
  

  cout << "Map file data read in." << endl;
}

//import gamma tracking waveform basis from a ROOT file
void GT_import_basis(TFile *coarse_basis_file, TFile *fine_basis_file, GT_basis *gt_basis){

  if(FINE_BASIS_BINFACTOR<COARSE_BASIS_BINFACTOR){
    cout << "ERROR: GT_import_basis - fine basis binning must be equal to or finer than the coarse basis binning." << endl;
    exit(-1);
  }
  if(fmod((FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR),1) > 0.){
    cout << "ERROR: GT_import_basis - the ratio FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR must be an integer." << endl;
    exit(-1);
  }

  //read in histograms from basis file
  if(!coarse_basis_file->IsOpen()) {
    cout << "ERROR: GT_import_basis - could not open coarse basis file!" << endl;
    exit(-1);
  }
  if(!fine_basis_file->IsOpen()) {
    cout << "ERROR: GT_import_basis - could not open fine basis file!" << endl;
    exit(-1);
  }

  cout << "Reading in coarse waveform basis data..." << endl;
  char hname[64];

  //setup histograms for the coarse basis
  if((gt_basis->basisHPCoarse = (TH1I*)coarse_basis_file->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the coarse waveform basis." << endl;
    exit(-1);
  }
  Int_t coarseBasisBinsR = VOXEL_BINS_R*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsZ = VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR;
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < coarseBasisBinsR; k++){
      Int_t numAngleBinsAtR = 4*getNumAngleBins(k,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
      for(int j = 0; j < numAngleBinsAtR; j++){
        for(int i = 0; i < coarseBasisBinsZ; i++){
          Int_t basisInd = l*coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ + k*coarseBasisBinsAngle*coarseBasisBinsZ + j*coarseBasisBinsZ + i;
          sprintf(hname,"basisPos%i_r%ito%i_angle%ito%i_z%ito%i",l,k*MAX_VAL_R/coarseBasisBinsR,(k+1)*MAX_VAL_R/coarseBasisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/coarseBasisBinsZ,(i+1)*MAX_VAL_Z/coarseBasisBinsZ);
          if((gt_basis->coarseBasis[basisInd] = (TH1D*)coarse_basis_file->Get(hname))==NULL){
            //cout << "No coarse waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
          }
          //list->Add(gt_basis->coarseBasis[basisInd]);
        }
      }
    }
  }
  
  cout << "Reading in fine waveform basis data..." << endl;
  //setup histograms for the basis 
  if((gt_basis->basisHPFine = (TH1I*)fine_basis_file->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the fine waveform basis." << endl;
    exit(-1);
  }
  Int_t fineBasisBinsR = VOXEL_BINS_R*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsZ = VOXEL_BINS_Z*FINE_BASIS_BINFACTOR;
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < fineBasisBinsR; k++){
      Int_t numAngleBinsAtR = 4*getNumAngleBins(k,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR; //x4 since covering 2pi rather than pi/2
      for(int j = 0; j < numAngleBinsAtR; j++){
        for(int i = 0; i < fineBasisBinsZ; i++){
          Int_t basisInd = l*fineBasisBinsR*fineBasisBinsAngle*fineBasisBinsZ + k*fineBasisBinsAngle*fineBasisBinsZ + j*fineBasisBinsZ + i;
          sprintf(hname,"basisPos%i_r%ito%i_angle%ito%i_z%ito%i",l,k*MAX_VAL_R/fineBasisBinsR,(k+1)*MAX_VAL_R/fineBasisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/fineBasisBinsZ,(i+1)*MAX_VAL_Z/fineBasisBinsZ);
          if((gt_basis->fineBasis[basisInd] = (TH1D*)fine_basis_file->Get(hname))==NULL){
            //cout << "No fine waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
          }
          //list->Add(gt_basis->fineBasis[basisInd]);
        }
      }
    }
  }

  cout << "All waveform basis data read in." << endl;
}

//rotate the sub-segment hit position vector into a position withina clover
TVector3 GT_transform_position_to_clover(TTigressHit *tigress_hit, TVector3 *hit_pos){

  TVector3 aPos(*hit_pos);

  //rotate hit position depending on core # (R,G,B,W)
  switch(tigress_hit->GetArrayNumber() % 4){
    case 3:
      //white
      aPos.RotateZ(TMath::Pi()/2.);
      aPos.SetX(aPos.X() - 27.5);
      aPos.SetY(aPos.Y() + 27.5);
      break;
    case 2:
      //red
      aPos.RotateZ(TMath::Pi());
      aPos.SetX(aPos.X() - 27.5);
      aPos.SetY(aPos.Y() - 27.5);
      break;
    case 1:
      //green
      aPos.RotateZ(3.*TMath::Pi()/2.);
      aPos.SetX(aPos.X() + 27.5);
      aPos.SetY(aPos.Y() - 27.5);
      break;
    case 0:
    default:
      //blue
      aPos.SetX(aPos.X() + 27.5);
      aPos.SetY(aPos.Y() + 27.5);
      break;
  }

  return aPos;
}

//transforms a hit position within a core to an absolute hit position
TVector3 GT_transform_position_to_absolute(TTigressHit *tigress_hit, TVector3 *hit_pos){

  TVector3 aPos(*hit_pos);

  //rotate hit position depending on core # (R,G,B,W)
  switch(tigress_hit->GetArrayNumber() % 4){
    case 3:
      //white
      aPos.RotateZ(TMath::Pi()/2.);
      break;
    case 2:
      //red
      aPos.RotateZ(TMath::Pi());
      break;
    case 1:
      //green
      aPos.RotateZ(3.*TMath::Pi()/2.);
      break;
    case 0:
    default:
      //blue
      break;
  }

  //express sub-segment hit position vector relative to the centre of the hit segment
  if(aPos.Z() > 30.){
    //segs 4-7
    aPos.SetZ(aPos.Z() - 45.);
  }else{
    //segs 0-3
    aPos.SetZ(aPos.Z() - 15.);
  }
  aPos.SetPerp(aPos.Perp() - 15.);

  //rotate the sub-segment hit position vector into the lab frame
  TVector3 direction = tigress_hit->GetPosition().Unit();
  aPos.RotateUz(direction); // direction must be TVector3 of unit length

  return tigress_hit->GetPosition() + aPos;
}

//returns the Doppler corrected energy from the sub-segment position
//source_beta: % speed of light that the source (beam/recoil) is travelling at (0-1)
//source_vec: vector for the source (beam/recoil) direction
//tigress_hit: GRSIsort TTigressHit object for the TIGRESS hit information
//hit_pos: the sub-segment hit position, from GT_get_pos_direct, GT_get_pos_gridsearch, or similar
double GT_get_doppler(double source_beta, TVector3 *source_vec = nullptr, TTigressHit *tigress_hit = nullptr, TVector3* hit_pos = nullptr){
  if(source_vec == nullptr){
    source_vec = tigress_hit->GetBeamDirection();
  }
  if(tigress_hit == nullptr){
    cout << "GT_get_doppler: TIGRESS hit is NULL pointer." << endl;
    return 0.;
  }else if(hit_pos == nullptr){
    cout << "GT_get_doppler: hit position is NULL pointer." << endl;
    return 0.;
  }
  double tmp   = 0;
  double gamma = 1 / (sqrt(1 - pow(source_beta, 2)));
  tmp = tigress_hit->GetEnergy() * gamma * (1 - source_beta * TMath::Cos(GT_transform_position_to_absolute(tigress_hit,hit_pos).Angle(*source_vec)));
  return tmp;
}

//get the maximum charge segment, returns -1 if no segment with max charge
Int_t getMaxChargeSegHit(TTigressHit *tigress_hit){

  const Double_t coreCharge = tigress_hit->GetCharge();
  Double_t maxSegCharge = 0.;
  Int_t maxChargeSeg = -1;
  Int_t segsInData = 0;
  Int_t numSamples = -1;
  bool goodWaveforms = true;
  if(tigress_hit->GetSegmentMultiplicity() == NSEG){
    for(int i = 0; i < NSEG; i++){
      if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.2*coreCharge){
        if(tigress_hit->GetSegmentHit(i).GetCharge() > maxSegCharge){
          maxSegCharge = tigress_hit->GetSegmentHit(i).GetCharge();
          maxChargeSeg = i;
        }
      }
    }
  }

  return maxChargeSeg;

}

//Get hit position using the direct method
//This method requires the tigress hit (which contains waveform information for all segments),
//as well as the index of the segment hit to map
TVector3 GT_get_pos_direct(TTigressHit *tigress_hit, GT_map *gt_map){

  //find the segment with the largest energy deposit, this will be the 'hit' segment
  //and check the number of segments with significant energy deposits, since this method only works when
  //one segment is hit
  //and check that the waveforms are the same size, and that there are no duplicate segments
  const Double_t coreCharge = tigress_hit->GetCharge();
  Double_t maxSegCharge = 0.;
  Int_t maxChargeSeg = -1;
  Int_t numSegHits = 0; //counter for the number of segments with a hit (ie. over the threshold energy)
  Int_t segsInData = 0;
  Int_t numSamples = -1;
  bool goodWaveforms = true;
  if(tigress_hit->GetSegmentMultiplicity() == NSEG){
    for(int i = 0; i < NSEG; i++){
      if(numSamples < 0){
        numSamples = tigress_hit->GetSegmentHit(i).GetWaveform()->size();
      }
      //make sure all segments in the data are different
      if(segsInData&(1<<(tigress_hit->GetSegmentHit(i).GetSegment()-1))){
        //cout << "Entry " << jentry << ", multiple hits in one segment." << endl;
        goodWaveforms = false;
        break;
      }else{
        segsInData|=(1<<(tigress_hit->GetSegmentHit(i).GetSegment()-1));
      }
      if(tigress_hit->GetSegmentHit(i).GetWaveform()->size()!=numSamples){
        //cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
        goodWaveforms = false;
        break;
      }
      if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.2*coreCharge){
        if(tigress_hit->GetSegmentHit(i).GetCharge() > maxSegCharge){
          maxSegCharge = tigress_hit->GetSegmentHit(i).GetCharge();
          maxChargeSeg = i;
        }
        numSegHits++;
      }
    }
  }

  if((goodWaveforms)&&(numSegHits==1)&&(maxChargeSeg >= 0)){

    Int_t arrayPos = tigress_hit->GetArrayNumber();
    if((arrayPos>=0)&&(arrayPos<(NPOS*NCORE))){

      //calculate all ordering parameters (see ordering_parameter_calc.cxx)
      double rho = calc_ordering(tigress_hit,maxChargeSeg,0);
      if(rho != BAD_RETURN){
        double phi = calc_ordering(tigress_hit,maxChargeSeg,1);
        if(phi != BAD_RETURN){
          double zeta = calc_ordering(tigress_hit,maxChargeSeg,2);
          if(zeta != BAD_RETURN){

            Int_t segNum = tigress_hit->GetSegmentHit(maxChargeSeg).GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

            //here is where the mapping happens
            double r=-1.;
            double angle=-1.;
            double z=-1.;
            if(gt_map->rMap[arrayPos*NSEG + segNum]!=NULL){
              r = gt_map->rMap[arrayPos*NSEG + segNum]->GetBinContent(gt_map->rMap[arrayPos*NSEG + segNum]->FindBin(rho));
            }
            Int_t rInd = (Int_t)(r*VOXEL_BINS_R/MAX_VAL_R);
            if(rInd < VOXEL_BINS_R){
              if(gt_map->angleMap[arrayPos*NSEG*VOXEL_BINS_R + segNum*VOXEL_BINS_R + rInd]!=NULL){
                angle = gt_map->angleMap[arrayPos*NSEG*VOXEL_BINS_R + segNum*VOXEL_BINS_R + rInd]->GetBinContent(gt_map->angleMap[arrayPos*NSEG*VOXEL_BINS_R + segNum*VOXEL_BINS_R + rInd]->FindBin(phi));
                if(angle>=MAX_VAL_ANGLE){
                  angle=MAX_VAL_ANGLE-0.001; //have seen rare events where angle is exactly 90 degrees
                }
                Int_t angleInd = (Int_t)(angle*getNumAngleBins(rInd,1.0)/(MAX_VAL_ANGLE*1.0));
                //cout << "seg: " << segNum << ", r index: " << rInd << ", angle: " << angle << ", angle index: " << angleInd << ", num angle bins: " << getNumAngleBins(rInd,1.0) << endl;
                if(gt_map->zMap[arrayPos*NSEG*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]!=NULL){
                  z = gt_map->zMap[arrayPos*NSEG*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]->GetBinContent(gt_map->zMap[arrayPos*NSEG*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]->FindBin(zeta));
                  if(z>=MAX_VAL_Z){
                    z=MAX_VAL_Z-0.001; //have seen rare events where z is exactly 90 mm
                  }

                  if((r>0.)&&(z>0.)&&(angle>0.)){
                    if(segNum<=3){
                      //r corresponds to the distance from the central contact at z=30
                      r = sqrt(r*r - (30.-z)*(30.-z));
                    }
                    if((r==r)&&(angle==angle)&&(z==z)){
                      angle += 90.0*(segNum%4);
                      //cout << "seg: " << segNum << ", angle: " << angle << endl;
                      return TVector3(r*cos(angle*M_PI/180.),r*sin(angle*M_PI/180.),z);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  //only get here if something failed
  return TVector3(BAD_RETURN,0,0);
}

//Get hit position using the adaptive grid search method
//This method requires the tigress hit (which contains waveform information for all segments)
//For hits in multiple segments, the hit position is computed in the segment with the largest energy deposit 
TVector3 GT_get_pos_gridsearch(TTigressHit *tigress_hit, GT_basis *gt_basis){

  if(randGen==NULL){
    randGen = new TRandom3();
  }

  const Double_t coreCharge = tigress_hit->GetCharge();
  const Int_t coarseBasisBinsR = VOXEL_BINS_R*COARSE_BASIS_BINFACTOR;
  const Int_t coarseBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR;
  const Int_t coarseBasisBinsZ = VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR;
  const Int_t fineBasisBinsR = VOXEL_BINS_R*FINE_BASIS_BINFACTOR;
  const Int_t fineBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR;
  const Int_t fineBasisBinsZ = VOXEL_BINS_Z*FINE_BASIS_BINFACTOR;
  const Int_t basisBinRatio = (Int_t)(FINE_BASIS_BINFACTOR/(COARSE_BASIS_BINFACTOR*1.0));
  const std::vector<Short_t> *segwf;

  if((coreCharge <= BASIS_MIN_ENERGY)||(coreCharge <= 0)) return TVector3(BAD_RETURN,0,0); //bad energy
  bool isHit = false;
  //cout << "Number of segments: " << tigress_hit->GetSegmentMultiplicity() << endl;
  Int_t numSegHits = 0; //counter for the number of segments with a hit (ie. over the threshold energy)
  Int_t evtSegHP = 0; //event segment hitpattern, which will be compared against hitpatterns in the basis
  if(tigress_hit->GetSegmentMultiplicity() == NSEG){
    //all segments have waveforms
    evtSegHP = 0;
    //check that the waveforms are the same size, and that there are no duplicate segments
    bool goodWaveforms = true;
    Int_t segsInData = 0;
    Double_t maxSegCharge = 0.;
    Int_t maxChargeSeg = -1;
    Int_t numSamples = -1;
    
    for(int i = 0; i < NSEG; i++){
      //make sure all segments in the data are different
      if(segsInData&(1<<(tigress_hit->GetSegmentHit(i).GetSegment()-1))){
        //cout << "Entry " << jentry << ", multiple hits in one segment." << endl;
        goodWaveforms = false;
        break;
      }else{
        segsInData|=(1<<(tigress_hit->GetSegmentHit(i).GetSegment()-1));
      }
      if(numSamples < 0){
        numSamples = tigress_hit->GetSegmentHit(i).GetWaveform()->size();
      }
      if(tigress_hit->GetSegmentHit(i).GetWaveform()->size()!=numSamples){
        //cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
        goodWaveforms = false;
        break;
      }
      if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.2*coreCharge){
        if(tigress_hit->GetSegmentHit(i).GetCharge() > maxSegCharge){
          maxSegCharge = tigress_hit->GetSegmentHit(i).GetCharge();
          maxChargeSeg = tigress_hit->GetSegmentHit(i).GetSegment()-1; //0-indexed
        }
        evtSegHP|=(1<<(tigress_hit->GetSegmentHit(i).GetSegment()-1)); //GetSegment() is 1-indexed
        numSegHits++;
      }
    }
    if(evtSegHP==0){
      goodWaveforms = false;
    }

    if(goodWaveforms){
      //cout << "maxChargeSeg: " << maxChargeSeg << endl;

      Int_t arrayPos = tigress_hit->GetArrayNumber();
      if((arrayPos>=0)&&(arrayPos<(NPOS*NCORE))){

        Double_t minChisqCoarse = BIG_NUMBER;
        Int_t minIndCoarse = -1;

        if(numSegHits<=NSEG){
          Int_t hitInd[NSEG];
          for(int i=0;i<NSEG;i++){
            hitInd[i]=-1;
          }
          for(int i=arrayPos*coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ; i<(arrayPos+1)*coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ; i++){
            if(gt_basis->coarseBasis[i]!=NULL){
              //scaled sum of multiple basis waveforms case
              //cout << "entry " << jentry << ", scaled sum case." << endl;
              Int_t hitSeg[NSEG];
              Double_t scaleFacHit[NSEG];
              hitSeg[0] = commonHitInHP(evtSegHP,(Int_t)gt_basis->basisHPCoarse->GetBinContent(i+1),NSEG);
              if(hitSeg[0]==maxChargeSeg){
                bool hitsFound=true;
                hitInd[0] = i;
                for(int hitNum=1; hitNum<numSegHits; hitNum++){
                  for(int i2=i+1; i2<(arrayPos+1)*coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ; i2++){
                    if(gt_basis->coarseBasis[i2]!=NULL){
                      hitSeg[hitNum] = commonHitInHP(evtSegHP,(Int_t)gt_basis->basisHPCoarse->GetBinContent(i2+1),NSEG);
                      if(hitSeg[hitNum]>=0){
                        bool goodSeg=true;
                        for(int k=0;k<hitNum;k++){
                          if(hitSeg[hitNum]==hitSeg[k]){
                            goodSeg=false;
                          }
                        }
                        if(goodSeg){
                          hitInd[hitNum] = i2;
                          break;
                        }
                      }
                    }
                  }
                  if(hitInd[hitNum] < 0){
                    //basis indices not found
                    hitsFound=false;
                    break;
                  }
                }
                
                if(hitsFound){
                  
                  //get scaling factors for each hit
                  for(int j=0; j<NSEG; j++){
                    Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment()-1; //0-indexed
                    for(int hitNum=0; hitNum<numSegHits; hitNum++){
                      if(segNum==hitSeg[hitNum]){
                        scaleFacHit[hitNum] = tigress_hit->GetSegmentHit(j).GetCharge()/coreCharge;
                      }
                    }
                    //cout << "Entry " << jentry << ", num seg hits: " << numSegHits << ", scaleFactor 0: " << scaleFacHit[0] << endl;
                  }

                  //compute chisq
                  Double_t chisq = 0.;
                  for(int j=0; j<NSEG; j++){
                    Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                    Double_t hitWeight = 1.0;
                    if(evtSegHP&(1<<(segNum-1))){
                      hitWeight = GRID_HIT_SEG_WEIGHT;
                    }else{
                      hitWeight = GRID_NONHIT_SEG_WEIGHT;
                    }
                    segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                    Double_t seg_waveform_baseline = 0.;
                    for(int k = 0; k < BASELINE_SAMPLES; k++){
                      seg_waveform_baseline += segwf->at(k);
                    }
                    seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                    for(int k=0; k<numSamples; k++){
                      Double_t wfrmSampleVal = fabs((segwf->at(k) - seg_waveform_baseline)/coreCharge);
                      Double_t basisSampleVal = 0.;
                      for(int hitNum=0; hitNum<numSegHits; hitNum++){
                        basisSampleVal += scaleFacHit[hitNum]*gt_basis->coarseBasis[hitInd[hitNum]]->GetBinContent(segNum*SAMPLES + k + 1);
                      }
                      chisq += pow(wfrmSampleVal - basisSampleVal,2)*hitWeight;
                    }
                  }

                  //check if chisq is at minimum
                  if((chisq>0)&&(chisq<minChisqCoarse)){
                    minChisqCoarse = chisq;
                    minIndCoarse = hitInd[0]; //save the index of the hit in the segment with the largest energy deposit
                  }

                }

              }
            }
          }
          //cout << "Entry " << jentry << " chisq val: " << minChisqCoarse << ", min indices: " << minIndCoarse[0] << ", " << minIndCoarse[1] << endl;
        }else{
          cout << "WARNING: grid search has unhandled number of hit segments (" << numSegHits << ")." << endl;
          return TVector3(BAD_RETURN,0,0);
        }

        //check that the coarse basis sort was successful
        if(minIndCoarse < 0){
          //at least one hit didn't have a basis index found
          return TVector3(BAD_RETURN,0,0);
        }


        //now do everything with the fine basis
        if(minChisqCoarse<MAX_BASIS_SORT_CHISQ){

          Int_t minRBinFine = -1;
          Int_t minAngleBinFine = -1;
          Int_t minZBinFine = -1;
          Int_t rBinFine, angleBinFine, zBinFine, basisInd;
          
          Double_t minChisqFine = BIG_NUMBER;

          //a best-fit coarse basis voxel was found, figure out the corresponding bin numbers
          Int_t rBinCoarse = (Int_t)(minIndCoarse/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
          Int_t angleBinCoarse = (Int_t)((minIndCoarse - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
          Int_t zBinCoarse = (minIndCoarse - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBinCoarse*coarseBasisBinsZ);
          //cout << "rBinCoarse: " << rBinCoarse << ", angleBinCoarse: " << angleBinCoarse << ", zBinCoarse: " << zBinCoarse << endl;

          for(int i=0; i<(basisBinRatio*basisBinRatio*basisBinRatio); i++){

            rBinFine = (rBinCoarse*basisBinRatio) + (Int_t)(i/(basisBinRatio*basisBinRatio*1.0));
            angleBinFine = (angleBinCoarse*basisBinRatio) + (((Int_t)(i/(basisBinRatio*1.0)) % basisBinRatio) );
            zBinFine = (zBinCoarse*basisBinRatio) + (i % basisBinRatio);
            basisInd = arrayPos*fineBasisBinsR*fineBasisBinsAngle*fineBasisBinsZ + rBinFine*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine*fineBasisBinsZ + zBinFine;

            if(gt_basis->fineBasis[basisInd]!=NULL){
            
              //signal resprented by either a single basis waveform with both segments hit (typically near segment boundaries)
              //or by multiple basis waveforms with hits in different segments
              //in the latter case, should scale each basis waveform based on the corresponding segment energy and then 
              //add together before performing chisq test
            
              //scaled sum of multiple basis waveforms case
              Int_t hitSeg;
              Double_t scaleFacHit;
              hitSeg = commonHitInHP(evtSegHP,(Int_t)gt_basis->basisHPFine->GetBinContent(basisInd+1),NSEG);
              if(hitSeg==maxChargeSeg){

                //get scaling factors for each hit
                for(int j=0; j<NSEG; j++){
                  Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment()-1; //0-indexed
                  if(segNum==hitSeg){
                    scaleFacHit = tigress_hit->GetSegmentHit(j).GetCharge()/coreCharge;
                    break;
                  }
                }

                //compute chisq
                Double_t chisq = 0.;
                for(int j=0; j<NSEG; j++){
                  Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                  Double_t hitWeight = 1.0;
                  if(evtSegHP&(1<<(segNum-1))){
                    hitWeight = GRID_HIT_SEG_WEIGHT;
                  }else{
                    hitWeight = GRID_NONHIT_SEG_WEIGHT;
                  }
                  segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                  Double_t seg_waveform_baseline = 0.;
                  for(int k = 0; k < BASELINE_SAMPLES; k++){
                    seg_waveform_baseline += segwf->at(k);
                  }
                  seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                  for(int k=0; k<numSamples; k++){
                    Double_t wfrmSampleVal = fabs((segwf->at(k) - seg_waveform_baseline)/coreCharge);
                    Double_t basisSampleVal = scaleFacHit*gt_basis->fineBasis[basisInd]->GetBinContent(segNum*SAMPLES + k + 1);
                    chisq += pow(wfrmSampleVal - basisSampleVal,2)*hitWeight;
                  }
                }

                //check if chisq is at minimum
                if((chisq>0)&&(chisq<minChisqFine)){
                  minChisqFine = chisq;
                  for(int hitNum=0; hitNum<numSegHits; hitNum++){
                    minRBinFine = rBinFine;
                    minAngleBinFine = angleBinFine;
                    minZBinFine = zBinFine;
                  }
                }
              }
            }
          }

          //a best-fit fine basis waveform was found, figure out the corresponding position and return it
          if(minChisqFine<MAX_BASIS_SORT_CHISQ){
            double rVal, angleVal, zVal;
            Int_t numAngleBinsAtR = 4*getNumAngleBins(minRBinFine,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR; //x4 since covering 2pi rather than pi/2
            rVal = (minRBinFine+randGen->Uniform(1.0))*MAX_VAL_R/(1.0*fineBasisBinsR);
            angleVal = (minAngleBinFine+randGen->Uniform(1.0))*360./(1.0*numAngleBinsAtR);
            zVal = (minZBinFine+randGen->Uniform(1.0))*MAX_VAL_Z/(1.0*fineBasisBinsZ);
            //cout << "  vals:" << rVal << " " << angleVal << " " << zVal << endl;
            return TVector3(rVal*cos(angleVal*M_PI/180.),rVal*sin(angleVal*M_PI/180.),zVal);
          }

        }
      }else{
        cout << "WARNING: Bad array position (" << arrayPos << ")" << endl;
      }
      
    }
  }

  //only get here if something failed
  return TVector3(BAD_RETURN,0,0);
      
}

