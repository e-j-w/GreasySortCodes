//function for calculation of ordering parameters
//this enforces the single segment hit condition assumed by the mapping process, using SEGMENT_ENERGY_THRESHOLD
double calc_ordering(TTigressHit * tigress_hit, const Int_t i, const Int_t jentry, const Int_t waveform_t0, const Int_t parameterNum) {

  //lists of adjacent segments in the TIGRESS array (zero-indexed)
  Int_t phiAdjSeg1[8] = {3,0,1,2,7,4,5,6};
  Int_t phiAdjSeg2[8] = {1,2,3,0,5,6,7,4};
  Int_t    zAdjSeg[8] = {4,5,6,7,0,1,2,3};
  
  double dno = 0.; //placeholder for denominator
  Int_t one;
  Int_t offset = 0;

  TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(i);

  if(segment_hit.GetCharge() < SEGMENT_ENERGY_THRESHOLD){
    return BAD_RETURN;
  }
  if(segment_hit.GetCharge() > MAX_ENERGY_SINGLE_INTERACTION){
    return BAD_RETURN; //energy too high - event likely to result from multiple hits
  }

  Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed
  //cout << "segment " << segNum << " energy: " << segment_hit.GetEnergy() << endl;

  const std::vector<Short_t> *segwf;
  segwf = segment_hit.GetWaveform();
  if(segwf->size() != SAMPLES){
    cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
    return BAD_RETURN;
  }
  

  if(parameterNum==0){
    //construct rho, the ordering parameter for the radius
    //see Eq. 4 of NIM A 729 (2013) 198-206
    double sampleAvg = 0.;
    //double term2 = 0.;
    double rho = 0.;
    const std::vector<Short_t> *segwf2, *segwf3, *segwf4;
    bool found1 = false;
    bool found2 = false;
    bool found3 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if((j!=i)&&(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD)){
        return BAD_RETURN;
      }
      //cout << "energy: " << tigress_hit->GetSegmentHit(j).GetEnergy() << ", charge: " << tigress_hit->GetSegmentHit(j).GetCharge() << endl;
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
        found1=true;
        segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
        found2=true;
        segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
        found3=true;
        segwf4 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
    }
    if((!found1)||(!found2)||(!found3)){
      cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
      return BAD_RETURN;
    }else if((segwf2->size() != SAMPLES)||(segwf3->size() != SAMPLES)||(segwf4->size() != SAMPLES)){
      cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    for(int j=1;j<SAMPLES-1;j++){
      sampleAvg += (j)*(segwf->at(j+1) - segwf->at(j-1))/2.0;
      //term2 += segwf2->at(j) + segwf3->at(j) + segwf4->at(j);
      dno += (segwf->at(j+1) - segwf->at(j-1))/2.0;
    }
    sampleAvg /= dno;
    for(int j=1;j<SAMPLES-1;j++){
      rho += pow(j - sampleAvg,3.0)*(segwf->at(j+1) - segwf->at(j-1))/2.0;// - 0.0005*term2;
    }
    rho /= dno;
    if((dno==0.)||(rho!=rho)){
      cout << "Entry " << jentry << ", cannot compute rho parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    //cout << "rho: " << rho << endl;
    return rho;
  }else if(parameterNum==1){
    //contruct phi, the ordering parameter for the angle
    //see Eq. 3 of NIM A 729 (2013) 198-206
    double phi = 0.;
    const std::vector<Short_t> *segwf2, *segwf3;
    bool found1 = false;
    bool found2 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if((j!=i)&&(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD)){
        return BAD_RETURN;
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
        found1=true;
        segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
        found2=true;
        segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
    }
    if((!found1)||(!found2)){
      cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute phi parameter." << endl;
      return BAD_RETURN;
    }else if((segwf2->size() != SAMPLES)||(segwf3->size() != SAMPLES)){
      cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    for(int j=0;j<SAMPLES;j++){
      phi += segwf2->at(j)*segwf2->at(j) - segwf3->at(j)*segwf3->at(j);
      dno += segwf2->at(j)*segwf2->at(j) + segwf3->at(j)*segwf3->at(j);
    }
    phi /= dno;
    if((dno==0.)||(phi!=phi)){
      cout << "Entry " << jentry << ", cannot compute phi parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    /*double maxVall = -1E30;
    double minVall = 1E30;
    double maxValr = -1E30;
    double minValr = 1E30;
    for(int j=0;j<WAVEFORM_SAMPLING_WINDOW ;j++){
      if(segwf2->at(waveform_t0+j) > maxVall){
        maxVall = segwf2->at(waveform_t0+j);
      }
      if(segwf2->at(waveform_t0+j) < minVall){
        minVall = segwf2->at(waveform_t0+j);
      }
      if(segwf3->at(waveform_t0+j) > maxValr){
        maxValr = segwf3->at(waveform_t0+j);
      }
      if(segwf3->at(waveform_t0+j) < minValr){
        minValr = segwf3->at(waveform_t0+j);
      }
    }
    if((minVall == 1E30)||(maxVall == 0)||(minValr == 1E30)||(maxValr == 0)){
      cout << "Entry " << jentry << ", cannot find maximum or minimum values for phi." << endl;
      cout << "vals: " << minVall << " " << maxVall << " " << minValr << " " << maxValr << endl;
      return BAD_RETURN;
    }
    phi = ((maxVall - minVall) - (maxValr - minValr))/((maxVall - minVall) + (maxValr - minValr));
    if((phi!=phi)){
      cout << "Entry " << jentry << ", cannot compute phi parameter (NaN)." << endl;
      return BAD_RETURN;
    }*/
    return phi;
  }else{
    //contruct zeta, the ordering parameter for the z direction
    //see Eq. 2 of NIM A 729 (2013) 198-206 (modified here)
    double zeta = 0.;
    const std::vector<Short_t> *segwf2;
    bool found2 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if((j!=i)&&(tigress_hit->GetSegmentHit(j).GetCharge() >= SEGMENT_ENERGY_NOHIT_THRESHOLD)){
        return BAD_RETURN;
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
        found2=true;
        segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
    }
    if(!found2){
      cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
      return BAD_RETURN;
    }else if(segwf2->size() != SAMPLES){
      cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    double maxVal = -1E30;
    double minVal = 1E30;
    for(int j=1;j<SAMPLES-1;j++){
      if(segwf2->at(j) > maxVal){
        maxVal = segwf2->at(j);
      }
      if(segwf2->at(j) < minVal){
        minVal = segwf2->at(j);
      }
      dno += (segwf->at(j+1) - segwf->at(j-1))/2.0;
    }
    if((minVal == 1E30)||(maxVal == 0)||(maxVal == minVal)){
      cout << "Entry " << jentry << ", cannot find maximum or minimum values for zeta." << endl;
      return BAD_RETURN;
    }
    zeta = (maxVal - minVal)/dno;
    if((dno==0.)||(zeta!=zeta)){
      cout << "Entry " << jentry << ", cannot compute zeta parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    //cout << "max: " << maxVal << ", min: " << minVal << ", dno: " << dno << ", zeta: " << zeta << endl;
    if(segNum>3){
      //back segment, reverse sign to make zeta increase with z
      zeta *= -1.;
      //zeta += ZETA_MAX;
    }//else{
    //  zeta -= ZETA_MAX;
    //}
    //cout << "zeta: " << zeta << endl;
    return zeta;
  }

}
