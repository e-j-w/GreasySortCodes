//function for calculation of ordering parameters
double calc_ordering(TTigressHit * tigress_hit, Int_t i, Int_t jentry, Int_t waveform_t0, Int_t parameterNum) {

  //lists of adjacent segments in the TIGRESS array (zero-indexed)
  Int_t phiAdjSeg1[8] = {3,0,1,2,7,4,5,6};
  Int_t phiAdjSeg2[8] = {1,2,3,0,5,6,7,4};
  Int_t    zAdjSeg[8] = {4,5,6,7,0,1,2,3};

  Int_t samples = 100; //number of samples per waveform
  Int_t sampling_window = 10; //number of waveform samples used to construct ordering parameters

  const std::vector<Short_t> *segwf, *segwf2, *segwf3;
  bool found1, found2;
  double dno;
  Int_t one;
  Int_t offset = 0;

  TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(i);

  Int_t posNum = tigress_hit->GetDetector()-1;
  Int_t coreNum = tigress_hit->GetCrystal();
  Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

  //cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
  segwf = segment_hit.GetWaveform();
  if((posNum < 0)||(posNum > 15)){
    cout << "Entry " << jentry << ", invalid array position: " << posNum << endl;
    return BAD_RETURN;
  }
  if(segwf->size() != samples){
    cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
    return BAD_RETURN;
  }


  if(parameterNum==0){
    //construct rho, the ordering parameter for the radius
    //see Eq. 4 of NIM A 729 (2013) 198-206
    double sampleAvg = 0.;
    double rho = 0.;
    dno = 0.; //placeholder for denominator
    for(int j=0;j<sampling_window;j++){
      sampleAvg += waveform_t0+j*(segwf->at(waveform_t0+j+1) - segwf->at(waveform_t0+j-1))/2.0;
      dno += (segwf->at(waveform_t0+j+1) - segwf->at(waveform_t0+j-1))/2.0;
    }
    sampleAvg /= dno;
    for(int j=0;j<sampling_window;j++){
      rho += pow(waveform_t0+j - sampleAvg,3.0)*(segwf->at(waveform_t0+j+1) - segwf->at(waveform_t0+j-1))/2.0;
    }
    rho /= dno;
    if((dno==0.)||(rho!=rho)){
      cout << "Entry " << jentry << ", cannot compute rho parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    return rho;
  }else if(parameterNum==1){
    //contruct phi, the ordering parameter for the angle
    //see Eq. 3 of NIM A 729 (2013) 198-206
    double phi = 0.;
    found1 = false;
    found2 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg1[segNum]){
        found1=true;
        segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == phiAdjSeg2[segNum]){
        found2=true;
        segwf2 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
    }
    if((!found1)||(!found2)){
      cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute phi parameter." << endl;
      cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
      return BAD_RETURN;
    }else if((segwf->size() != samples)||(segwf2->size() != samples)){
      cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    for(int j=0;j<sampling_window;j++){
      phi += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
      dno += segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j);
    }
    phi /= dno;
    if((dno==0.)||(phi!=phi)){
      cout << "Entry " << jentry << ", cannot compute phi parameter (NaN)." << endl;
      return BAD_RETURN;
    }
    return phi;
  }else{
    //contruct zeta, the ordering parameter for the z direction
    //see Eq. 2 of NIM A 729 (2013) 198-206 (modified here)
    double zeta = 0.;
    found1 = false;
    for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
      if(tigress_hit->GetSegmentHit(j).GetSegment()-1 == zAdjSeg[segNum]){
        found1=true;
        segwf3 = tigress_hit->GetSegmentHit(j).GetWaveform();
      }
    }
    if(!found1){
      cout << "Entry " << jentry << ", cannot get neighbouring segment wavefoms to compute zeta parameter." << endl;
      cout << "Entry " << jentry << ", position: " << posNum << ", core: " << coreNum << ", segment: " << segNum << endl;
      return BAD_RETURN;
    }else if(segwf3->size() != samples){
      cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
      return BAD_RETURN;
    }
    for(int j=0;j<sampling_window;j++){
      zeta += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) - segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) - segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
      dno += 2.0*segwf3->at(waveform_t0+j)*segwf3->at(waveform_t0+j) + segwf2->at(waveform_t0+j)*segwf2->at(waveform_t0+j) + segwf->at(waveform_t0+j)*segwf->at(waveform_t0+j);
    }
    zeta /= dno;
    if((dno==0.)||(zeta!=zeta)){
      cout << "Entry " << jentry << ", cannot compute zeta parameter (NaN)." << endl;
      return BAD_RETURN;
    }

    if(segNum>3){
      //back segment, reverse sign to make zeta increase with z
      zeta *= -1.;
    }
    return zeta;
  }

}
