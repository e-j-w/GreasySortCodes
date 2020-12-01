#include "common.h" //define all global variables here!


TH1D *coarseBasis[(Int_t)(VOXEL_BINS_R*COARSE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR*VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR)];
TH1D *fineBasis[(Int_t)(VOXEL_BINS_R*FINE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR)];
TRandom3 *randGen;

void fillAtHitInd(const Int_t hitInd, const Int_t fineBasisBinsR, const Int_t fineBasisBinsAngle, const Int_t fineBasisBinsZ, TRandom3 *randGen, TH3D *pos3DMap){
  Int_t rBinFine, angleBinFine, zBinFine;
  double rVal, angleVal, zVal;
  if(hitInd >= 0){
    //a best-fit fine basis waveform was found, figure out the corresponding bin
    rBinFine = (Int_t)(hitInd/(fineBasisBinsAngle*fineBasisBinsZ*1.0));
    angleBinFine = (Int_t)((hitInd - rBinFine*(fineBasisBinsAngle*fineBasisBinsZ*1.0))/(fineBasisBinsZ*1.0));
    zBinFine = (hitInd - rBinFine*(fineBasisBinsAngle*fineBasisBinsZ) - angleBinFine*fineBasisBinsZ);
    Int_t numAngleBinsAtR = 4*getNumAngleBins(rBinFine,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR; //x4 since covering 2pi rather than pi/2
    rVal = (rBinFine+randGen->Uniform(1.0))*MAX_VAL_R/(1.0*fineBasisBinsR);
    angleVal = (angleBinFine+randGen->Uniform(1.0))*360./(1.0*numAngleBinsAtR);
    zVal = (zBinFine+randGen->Uniform(1.0))*MAX_VAL_Z/(1.0*fineBasisBinsZ);
    //cout << "Filling at: x=" << rVal*cos(angleVal*M_PI/180.) << ", y=" << rVal*sin(angleVal*M_PI/180.) << ", z=" << zVal << ", hitInd: " << hitInd << endl;
    pos3DMap->Fill(rVal*cos(angleVal*M_PI/180.),rVal*sin(angleVal*M_PI/180.),zVal);
  }
}
void fillXYAtHitInd(const Int_t hitInd, const Int_t fineBasisBinsR, const Int_t fineBasisBinsAngle, const Int_t fineBasisBinsZ, TH2D *posXYMap){
  Int_t rBinFine, angleBinFine;
  double rVal, angleVal;
  if(hitInd >= 0){
    //a best-fit fine basis waveform was found, figure out the corresponding bin
    rBinFine = (Int_t)(hitInd/(fineBasisBinsAngle*fineBasisBinsZ*1.0));
    angleBinFine = (Int_t)((hitInd - rBinFine*(fineBasisBinsAngle*fineBasisBinsZ*1.0))/(fineBasisBinsZ*1.0));
    Int_t numAngleBinsAtR = 4*getNumAngleBins(rBinFine,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR; //x4 since covering 2pi rather than pi/2
    rVal = (rBinFine+0.5)*MAX_VAL_R/(1.0*fineBasisBinsR);
    angleVal = (angleBinFine+0.5)*360./(1.0*numAngleBinsAtR);
    posXYMap->Fill(rVal*cos(angleVal*M_PI/180.),rVal*sin(angleVal*M_PI/180.));
  }
}
void fillAtHitIndBasis(const Int_t hitInd, const Int_t basisBinsR, const Int_t basisBinsAngle, const Int_t basisBinsZ, TH3D *basisInd3DMap){
  Int_t rBinFine, angleBinFine, zBinFine;
  if(hitInd >= 0){
    //a best-fit fine basis waveform was found, figure out the corresponding bin
    rBinFine = (Int_t)(hitInd/(basisBinsAngle*basisBinsZ*1.0));
    angleBinFine = (Int_t)((hitInd - rBinFine*(basisBinsAngle*basisBinsZ*1.0))/(basisBinsZ*1.0));
    zBinFine = (hitInd - rBinFine*(basisBinsAngle*basisBinsZ) - angleBinFine*basisBinsZ);
    basisInd3DMap->Fill(rBinFine,angleBinFine,zBinFine);
  }
}

//checks if there is a common hit in the hitpatterns,
//and returns the bit index of the first common hit
int32_t commonHitInHP(int32_t hp1, int32_t hp2, int32_t max_search){
  if(max_search>32){
    cout << "WARNING: trying to search past the bounds of a 32-bit integer!" << endl;
    return -1;
  }
  uint32_t one = 1;
  for(int i=0;i<max_search;i++){
    if((hp1&(one<<i))&&(hp2&(one<<i))){
      return i;
    }
  }

  return -1;
}

void sortData(TFile *inputfile, const char *calfile, TH1I *basisHPCoarse, TH1I *basisHPFine, TH3D *coarseBasisInd3DMap, TH3D *basisInd3DMap, TH1I *numSegHitsHist, TH3D *pos3DMap, TH2D *posXYMap, TH3D *pos3DMap1Hit, TH3D *pos3DMap2Hit, TH1D *coarseChisqHist, TH1D *fineChisqHist){
  
  Int_t coarseBasisBinsR = VOXEL_BINS_R*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsZ = VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR;
  Int_t fineBasisBinsR = VOXEL_BINS_R*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsZ = VOXEL_BINS_Z*FINE_BASIS_BINFACTOR;

  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  cout << AnalysisTree->GetNtrees() << " tree files, details:" << endl;
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  cout << "Reading calibration file: " << calfile << endl;
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TTigress * tigress = 0;
  TTigressHit * tigress_hit;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    exit(-1);
  }

  const std::vector<Short_t> *segwf;
  const Int_t one = 1;
  const Int_t basisBinRatio = (Int_t)(FINE_BASIS_BINFACTOR/(COARSE_BASIS_BINFACTOR*1.0));
  //for (int jentry = 0; jentry < 100000; jentry++) {
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;
      Double_t coreCharge = tigress_hit->GetCharge();
      if(coreCharge <= 0) continue; //bad energy
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
        for(int i = 0; i < NSEG; i++){
          //make sure all segments in the data are different
          if(segsInData&(one<<(tigress_hit->GetSegmentHit(i).GetSegment()-1))){
            cout << "Entry " << jentry << ", multiple hits in one segment." << endl;
            goodWaveforms = false;
            break;
          }else{
            segsInData|=(one<<(tigress_hit->GetSegmentHit(i).GetSegment()-1));
          }
          if(tigress_hit->GetSegmentHit(i).GetWaveform()->size()!=SAMPLES){
            cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
            goodWaveforms = false;
            break;
          }
          if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.2*coreCharge){
            evtSegHP|=(one<<(tigress_hit->GetSegmentHit(i).GetSegment()-1)); //GetSegment() is 1-indexed
            numSegHits++;
          }
        }
        //cout << "HP: " << evtSegHP << ", num seg hits: " << numSegHits << endl;
        //cout << "good: " << goodWaveforms << ", num seg hits: " << numSegHits << endl;
        if(goodWaveforms){

          Double_t minChisqCoarse = BIG_NUMBER;
          
          Int_t minIndCoarse[2];
          minIndCoarse[0] = -1;
          minIndCoarse[1] = -1;

          //take action depending on the number of hit segments
          if(numSegHits==1){
            //signal best represented by a basis waveform from a position within the segment of interest
            for(int i=0; i<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i++){
              if(coarseBasis[i]!=NULL){
                //check that the event hitpattern matches the hitpattern in the basis
                //if(commonHitInHP(evtSegHP,(Int_t)basisHPCoarse->GetBinContent(i+1),8)>=0){
                if((Int_t)basisHPCoarse->GetBinContent(i+1) == evtSegHP){
                  //cout << "entry " << jentry << ", index " << i << ", single hit" << endl;
                  Double_t chisq = 0;
                  for(int j=0; j<NSEG; j++){
                    Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                    Double_t hitWeight = 1.0;
                    if(evtSegHP&(one<<(segNum-1))){
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
                    for(int k=0; k<SAMPLES; k++){
                      Double_t wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/coreCharge;
                      Double_t basisSampleVal = coarseBasis[i]->GetBinContent(segNum*SAMPLES + k + 1);
                      chisq += pow(wfrmSampleVal - basisSampleVal,2)*hitWeight;
                      /*if(evtSegHP==2){
                        rBin = (Int_t)(i/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                        angleBin = (Int_t)((i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                        zBin = (i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBin*coarseBasisBinsZ);
                        cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                        cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                        cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                      }*/
                    }
                  }
                  if((chisq>0.)&&(chisq<minChisqCoarse)){
                    minChisqCoarse = chisq;
                    minIndCoarse[0] = i;
                    minIndCoarse[1] = -1; //no 2nd basis waveform
                    //cout << "entry " << jentry << ", hit " << hitInd << ", basis index " << i << ", basis HP: " << basisHPCoarse->GetBinContent(i+1) << ", event HP: " << evtSegHP << ", chisq: " << chisq << endl;
                  }
                }
              }
            }
          }else if(numSegHits==2){
            //signal resprented by either a single basis waveform with both segments hit (typically near segment boundaries)
            //or by two basis waveforms with hits in different segments
            //in the latter case, should scale each basis waveform based on the corresponding segment energy and then 
            //add together before performing chisq test
            continue;
            
            for(int i=0; i<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i++){
              if(coarseBasis[i]!=NULL){
                //identical hitpattern case
                if((Int_t)basisHPCoarse->GetBinContent(i+1) == evtSegHP){
                  Double_t chisq = 0;
                  for(int j=0; j<NSEG; j++){
                    Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                    segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                    Double_t seg_waveform_baseline = 0.;
                    for(int k = 0; k < BASELINE_SAMPLES; k++){
                      seg_waveform_baseline += segwf->at(k);
                    }
                    seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                    for(int k=0; k<SAMPLES; k++){
                      Double_t wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/coreCharge;
                      Double_t basisSampleVal = coarseBasis[i]->GetBinContent(segNum*SAMPLES + k + 1);
                      chisq += pow(wfrmSampleVal - basisSampleVal,2);
                      /*if(evtSegHP==2){
                        rBin = (Int_t)(i/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                        angleBin = (Int_t)((i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                        zBin = (i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBin*coarseBasisBinsZ);
                        cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                        cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                        cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                      }*/
                    }
                  }
                  if((chisq>0)&&(chisq<minChisqCoarse)){
                    minChisqCoarse = chisq;
                    minIndCoarse[0] = i;
                    minIndCoarse[1] = -1; //no 2nd basis waveform
                    //cout << "index " << i << ", basis HP: " << basisHPCoarse->GetBinContent(i+1) << ", event HP: " << evtSegHP << ", chisq: " << chisq << endl;
                  }
                }else{
                  //scaled sum of 2 basis waveforms case
                  Int_t hitSeg1, hitSeg2;
                  Double_t scaleFacHit1, scaleFacHit2;
                  hitSeg1 = commonHitInHP(evtSegHP,(Int_t)basisHPCoarse->GetBinContent(i+1),8);
                  if(hitSeg1>=0){
                    for(int i2=i+1; i2<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i2++){
                      if(coarseBasis[i2]!=NULL){
                        hitSeg2 = commonHitInHP(evtSegHP,(Int_t)basisHPCoarse->GetBinContent(i2+1),8);
                        if((hitSeg2>=0)&&(hitSeg2!=hitSeg1)){
                          //we have two basis indices which are different and which summed together contain hits
                          //in the two segments seen in the real event
                          //cout << "entry: " << jentry << ", evtHP: " << evtSegHP << ", hit segments: [" << hitSeg1 << " " << hitSeg2 << "]" << endl;
                          Double_t chisq = 0;
                          for(int j=0; j<NSEG; j++){
                            Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment()-1; //0-indexed
                            if(tigress_hit->GetSegmentHit(j).GetCharge()>0.){
                              if(segNum==hitSeg1){
                                scaleFacHit1=coreCharge/tigress_hit->GetSegmentHit(j).GetCharge();
                              }else if(segNum==hitSeg2){
                                scaleFacHit2=coreCharge/tigress_hit->GetSegmentHit(j).GetCharge();
                              }
                            }
                          }
                          if((scaleFacHit1!=0.)&&(scaleFacHit2!=0.)){
                            for(int j=0; j<NSEG; j++){
                              Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                              segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                              Double_t seg_waveform_baseline = 0.;
                              for(int k = 0; k < BASELINE_SAMPLES; k++){
                                seg_waveform_baseline += segwf->at(k);
                              }
                              seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                              for(int k=0; k<SAMPLES; k++){
                                Double_t wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/coreCharge;
                                Double_t basisSampleVal = scaleFacHit1*coarseBasis[i]->GetBinContent(segNum*SAMPLES + k + 1);
                                Double_t basisSampleVal2 = scaleFacHit2*coarseBasis[i2]->GetBinContent(segNum*SAMPLES + k + 1);
                                chisq += pow(wfrmSampleVal - basisSampleVal - basisSampleVal2,2);
                                /*if(evtSegHP==2){
                                  rBin = (Int_t)(i/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                                  angleBin = (Int_t)((i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                                  zBin = (i - rBin*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBin*coarseBasisBinsZ);
                                  cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                                  cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                                  cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                                }*/
                              }
                              
                            }
                            if((chisq>0)&&(chisq<MAX_BASIS_SORT_CHISQ)&&(chisq<minChisqCoarse)){
                              minChisqCoarse = chisq;
                              minIndCoarse[0] = i;
                              minIndCoarse[1] = i2;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }else{
            //cout << "Entry " << jentry << " contains an unhandled number of hit segments (" << numSegHits << ")." << endl;
            continue;
          }

          coarseChisqHist->Fill(minChisqCoarse);

          //now do everything with the fine basis
          if(minChisqCoarse<MAX_BASIS_SORT_CHISQ){
            for(int cbh=0;cbh<2;cbh++){

              Double_t minChisqFine = BIG_NUMBER;

              Int_t minIndFine[2];
              minIndFine[0] = -1;
              minIndFine[1] = -1;

              //cout << "minIndCoarse: " << minIndCoarse[cbh] << endl;

              if(minIndCoarse[cbh] >= 0){

                fillAtHitIndBasis(minIndCoarse[cbh],coarseBasisBinsR,coarseBasisBinsAngle,coarseBasisBinsZ,coarseBasisInd3DMap);

                //a best-fit coarse basis waveform was found, figure out the corresponding bin
                Int_t rBinCoarse = (Int_t)(minIndCoarse[cbh]/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                Int_t angleBinCoarse = (Int_t)((minIndCoarse[cbh] - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                Int_t zBinCoarse = (minIndCoarse[cbh] - rBinCoarse*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBinCoarse*coarseBasisBinsZ);
                //rVal = (rBin+0.5)*MAX_VAL_R/(1.0*coarseBasisBinsR);
                //angleVal = (angleBin+0.5)*360./(1.0*coarseBasisBinsAngle);
                //zVal = (zBin+0.5)*MAX_VAL_Z/(1.0*coarseBasisBinsZ);
                /*if(rBinCoarse==0){
                  cout << "coarse basis bins: [ " << rBinCoarse << " " << angleBinCoarse << " " << zBinCoarse << " ]" << endl;
                  cout << "angle bins coarse: " << 4*getNumAngleBins(rBinCoarse,COARSE_BASIS_BINFACTOR) << endl;
                }*/
                
                

                for(int i=0; i<(basisBinRatio*basisBinRatio*basisBinRatio); i++){

                  Int_t rBinFine = (rBinCoarse*basisBinRatio) + (Int_t)(i/(basisBinRatio*basisBinRatio*1.0));
                  Int_t angleBinFine = (angleBinCoarse*basisBinRatio) + (((Int_t)(i/(basisBinRatio*1.0)) % basisBinRatio) );
                  Int_t zBinFine = (zBinCoarse*basisBinRatio) + (i % basisBinRatio);
                  Int_t basisInd = rBinFine*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine*fineBasisBinsZ + zBinFine;
                  //cout << "coarse basis bins [" << rBinCoarse << " " << angleBinCoarse << " " << zBinCoarse << "]" << endl;
                  //cout << "fine basis bin: " << basisInd << endl;
                  /*if(rBinCoarse==0){
                    cout << "fine basis bins [" << rBinFine << " " << angleBinFine << " " << zBinFine << "], angle bins fine : " << 4*getNumAngleBins(rBinFine,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR << ", fine basis index: " << basisInd << endl;
                  }*/
                  


                  if(fineBasis[basisInd]!=NULL){
                    //cout << "   fine basis index exists" << endl;
                    if(numSegHits>=1){
                      //signal may be represented by a basis waveform from a position within the segment of interest

                      //check that the event hitpattern matches the hitpattern in the basis
                      //cout << "Fine basis HP: " << basisHPFine->GetBinContent(basisInd+1) << ", event HP: " << evtSegHP << endl;
                      //if(commonHitInHP(evtSegHP,(Int_t)basisHPFine->GetBinContent(basisInd+1),8)>=0){
                      if(basisHPFine->GetBinContent(basisInd+1) == evtSegHP){
                        //cout << "Fine basis HP: " << basisHPFine->GetBinContent(basisInd+1) << endl;
                        Double_t chisq = 0;
                        for(int j=0; j<NSEG; j++){
                          Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                          Double_t hitWeight = 1.0;
                          if(evtSegHP&(one<<(segNum-1))){
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
                          for(int k=0; k<SAMPLES; k++){
                            Double_t wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/coreCharge;
                            Double_t basisSampleVal = fineBasis[basisInd]->GetBinContent(segNum*SAMPLES + k + 1);
                            chisq += pow(wfrmSampleVal - basisSampleVal,2)*hitWeight;
                            /*if(evtSegHP==2){
                              rBin = (Int_t)(i/(fineBasisBinsAngle*fineBasisBinsZ*1.0));
                              angleBin = (Int_t)((i - rBin*(fineBasisBinsAngle*fineBasisBinsZ*1.0))/(fineBasisBinsZ*1.0));
                              zBin = (i - rBin*(fineBasisBinsAngle*fineBasisBinsZ) - angleBin*fineBasisBinsZ);
                              cout << "HP: " << evtSegHP << ", rBin: " << rBin << ", angleBin: " << angleBin << ", zBin: " << zBin << ", segment: "; 
                              cout << tigress_hit->GetSegmentHit(j).GetSegment() << ", sample: " << k << ", sample val: ";
                              cout <<  wfrmSampleVal << ", basis sample val: " << basisSampleVal << endl;
                            }*/
                          }
                        }
                        //cout << "entry " << jentry <<  ", index " << basisInd << ", basis HP: " << basisHPFine->GetBinContent(basisInd+1) << ", event HP: " << evtSegHP << ", chisq: " << chisq << endl;
                        if((chisq>0)&&(chisq<minChisqFine)){
                          minChisqFine = chisq;
                          minIndFine[0] = basisInd;
                          minIndFine[1] = -1; //no 2nd basis waveform
                          //cout << "entry " << jentry <<  ", index " << basisInd << ", basis HP: " << basisHPFine->GetBinContent(basisInd+1) << ", event HP: " << evtSegHP << ", chisq: " << chisq << endl;
                          //cout << "min chisq" << endl;
                        }
                      }
                    }
                    /*if(numSegHits==2){
                      //signal resprented by either a single basis waveform with both segments hit (typically near segment boundaries)
                      //or by two basis waveforms with hits in different segments
                      //in the latter case, should scale each basis waveform based on the corresponding segment energy and then 
                      //add together before performing chisq test
                    
                      //scaled sum of 2 basis waveforms case
                      Int_t hitSeg1, hitSeg2;
                      Double_t scaleFacHit1, scaleFacHit2;
                      hitSeg1 = commonHitInHP(evtSegHP,(Int_t)basisHPFine->GetBinContent(basisInd+1),8);
                      if(hitSeg1>=0){
                        for(int i2=i+1; i2<(basisBinRatio*basisBinRatio*basisBinRatio); i2++){

                          rBinFine = (rBinCoarse*basisBinRatio) + (Int_t)(i2/(basisBinRatio*basisBinRatio*1.0));
                          angleBinFine = (angleBinCoarse*basisBinRatio) + (((Int_t)(i2/(basisBinRatio*1.0)) % basisBinRatio) );
                          zBinFine = (zBinCoarse*basisBinRatio) + (i2 % basisBinRatio);
                          Int_t basisInd2 = rBinFine*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine*fineBasisBinsZ + zBinFine;

                          if(fineBasis[basisInd2]!=NULL){
                            hitSeg2 = commonHitInHP(evtSegHP,(Int_t)basisHPFine->GetBinContent(basisInd2+1),8);
                            if((hitSeg2>=0)&&(hitSeg2!=hitSeg1)){
                              //we have two basis indices which are different and which summed together contain hits
                              //in the two segments seen in the real event
                              //cout << "entry: " << jentry << ", evtHP: " << evtSegHP << ", hit segments: [" << hitSeg1 << " " << hitSeg2 << "]" << endl;
                              Double_t chisq = 0;
                              for(int j=0; j<NSEG; j++){
                                Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment()-1; //0-indexed
                                if(tigress_hit->GetSegmentHit(j).GetCharge()>0.){
                                  if(segNum==hitSeg1){
                                    scaleFacHit1=coreCharge/tigress_hit->GetSegmentHit(j).GetCharge();
                                  }else if(segNum==hitSeg2){
                                    scaleFacHit2=coreCharge/tigress_hit->GetSegmentHit(j).GetCharge();
                                  }
                                }
                              }
                              if((scaleFacHit1!=0.)&&(scaleFacHit2!=0.)){
                                for(int j=0; j<NSEG; j++){
                                  Int_t segNum = tigress_hit->GetSegmentHit(j).GetSegment(); //1-indexed
                                  segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                                  Double_t seg_waveform_baseline = 0.;
                                  for(int k = 0; k < BASELINE_SAMPLES; k++){
                                    seg_waveform_baseline += segwf->at(k);
                                  }
                                  seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                                  for(int k=0; k<SAMPLES; k++){
                                    Double_t wfrmSampleVal = (segwf->at(k) - seg_waveform_baseline)/coreCharge;
                                    Double_t basisSampleVal = scaleFacHit1*fineBasis[basisInd]->GetBinContent(segNum*SAMPLES + k + 1);
                                    Double_t basisSampleVal2 = scaleFacHit2*fineBasis[basisInd2]->GetBinContent(segNum*SAMPLES + k + 1);
                                    chisq += pow(wfrmSampleVal - basisSampleVal - basisSampleVal2,2);
                                  }
                                }
                                if((chisq>0)&&(chisq<minChisqFine)){
                                  minChisqFine = chisq;
                                  minIndFine[0] = basisInd;
                                  minIndFine[1] = basisInd2;
                                }
                              }
                            }
                          }
                        }
                      }
                      
                    }*/
                  }
                }
              }

              fineChisqHist->Fill(minChisqFine);

              //fill histograms
              if(minChisqFine<MAX_BASIS_SORT_CHISQ){
                for(int i=0;i<2;i++){
                  if(minIndFine[i]>=0){
                    numSegHitsHist->Fill(numSegHits);
                    fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap);
                    fillAtHitIndBasis(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,basisInd3DMap);
                    fillXYAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,posXYMap);
                    if(numSegHits==1){
                      fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap1Hit);
                    }else if(numSegHits==2){
                      fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap2Hit);
                    }
                  }
                }
              }
            }
          }
          
        }
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
}

//Function which sorts hit positions using a pre-generated waveform basis
void sort_from_basis(const char *infile, const char *basisfileCoarse, const char *basisfileFine, const char *calfile, const char *outfile, bool inpList) {

  if(FINE_BASIS_BINFACTOR<COARSE_BASIS_BINFACTOR){
    cout << "ERROR: fine basis binning must be equal to or finer than the coarse basis binning." << endl;
    exit(-1);
  }
  if(fmod((FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR),1) > 0.){
    cout << "ERROR: the ratio FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR must be an integer." << endl;
    exit(-1);
  }

  randGen = new TRandom3();

  TList *list = new TList;

  //read in histograms from basis file
  TFile *inp = new TFile(basisfileCoarse,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open basis file!" << endl;
    exit(-1);
  }else{
    cout << "Opened waveform basis file: " << basisfileCoarse << endl;
  }
  
  //setup histograms for the basis 
  TH1I *basisHPCoarse, *basisHPFine;
  if((basisHPCoarse = (TH1I*)inp->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the coarse waveform basis." << endl;
    exit(-1);
  }
  Int_t coarseBasisBinsR = VOXEL_BINS_R*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR;
  Int_t coarseBasisBinsZ = VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR;
  for(int k = 0; k < coarseBasisBinsR; k++){
    Int_t numAngleBinsAtR = 4*getNumAngleBins(k,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
    for(int j = 0; j < numAngleBinsAtR; j++){
      for(int i = 0; i < coarseBasisBinsZ; i++){
        Int_t basisInd = k*coarseBasisBinsAngle*coarseBasisBinsZ + j*coarseBasisBinsZ + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/coarseBasisBinsR,(k+1)*MAX_VAL_R/coarseBasisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/coarseBasisBinsZ,(i+1)*MAX_VAL_Z/coarseBasisBinsZ);
        if((coarseBasis[basisInd] = (TH1D*)inp->Get(hname))==NULL){
          cout << "No coarse waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
          //cout << hname << endl;
          //getc(stdin);
        }
        //list->Add(coarseBasis[basisInd]);
      }
    }
  }
  cout << "Coarse waveform basis data read in." << endl;
  

  //read in histograms from basis file
  TFile *inp2 = new TFile(basisfileFine,"read");
  if (!inp2->IsOpen()) {
    cout << "ERROR: Could not open basis file!" << endl;
    exit(-1);
  }else{
    cout << "Opened waveform basis file: " << basisfileFine << endl;
  }
  
  //setup histograms for the basis 
  if((basisHPFine = (TH1I*)inp2->Get("basis_hitpattern"))==NULL){
    cout << "ERROR: no hitpattern in the fine waveform basis." << endl;
    exit(-1);
  }
  Int_t fineBasisBinsR = VOXEL_BINS_R*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR;
  Int_t fineBasisBinsZ = VOXEL_BINS_Z*FINE_BASIS_BINFACTOR;
  for(int k = 0; k < fineBasisBinsR; k++){
    Int_t numAngleBinsAtR = 4*getNumAngleBins(k,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR; //x4 since covering 2pi rather than pi/2
    for(int j = 0; j < numAngleBinsAtR; j++){
      for(int i = 0; i < fineBasisBinsZ; i++){
        Int_t basisInd = k*fineBasisBinsAngle*fineBasisBinsZ + j*fineBasisBinsZ + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/fineBasisBinsR,(k+1)*MAX_VAL_R/fineBasisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/fineBasisBinsZ,(i+1)*MAX_VAL_Z/fineBasisBinsZ);
        if((fineBasis[basisInd] = (TH1D*)inp2->Get(hname))==NULL){
          cout << "No fine waveform basis data for radial bin " << k << ", angle bin " << j << ", z bin " << i << endl;
        }
        //list->Add(fineBasis[basisInd]);
      }
    }
  }
  cout << "Fine waveform basis data read in." << endl;

  TH1I *numSegHitsHist = new TH1I("numSegHits","numSegHits",10,0,10);
  list->Add(numSegHitsHist);
  TH3D *pos3DMap = new TH3D("pos3DMap","pos3DMap",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);
  TH3D *pos3DMap1Hit = new TH3D("pos3DMap1Hit","pos3DMap1Hit",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap1Hit);
  TH3D *pos3DMap2Hit = new TH3D("pos3DMap2Hit","pos3DMap2Hit",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap2Hit);
  TH2D *posXYMap = new TH2D("posXYMap","posXYMap",40,-40,40,40,-40,40);
  list->Add(posXYMap);
  TH3D *coarseBasisInd3DMap = new TH3D("coarseBasisInd3DMap","coarseBasisInd3DMap",coarseBasisBinsR,0,coarseBasisBinsR,coarseBasisBinsAngle,0,coarseBasisBinsAngle,coarseBasisBinsZ,0,coarseBasisBinsZ);
  list->Add(coarseBasisInd3DMap);
  TH3D *basisInd3DMap = new TH3D("basisInd3DMap","basisInd3DMap",fineBasisBinsR,0,fineBasisBinsR,fineBasisBinsAngle,0,fineBasisBinsAngle,fineBasisBinsZ,0,fineBasisBinsZ);
  list->Add(basisInd3DMap);
  TH1D *coarseChisqHist = new TH1D("coarseChisqHist","coarseChisqHist",1000,0,10);
  list->Add(coarseChisqHist);
  TH1D *fineChisqHist = new TH1D("fineChisqHist","fineChisqHist",1000,0,10);
  list->Add(fineChisqHist);

  if(inpList){
    FILE *listFile;
    char name[256];
    if((listFile=fopen(infile,"r"))==NULL){
      cout << "ERROR: Could not open analysis tree list file!" << endl;
      exit(-1);
    }
    while(fscanf(listFile,"%s",name)!=EOF){
      TFile * inputfile = new TFile(name, "READ");
      if (!inputfile->IsOpen()) {
        cout << "ERROR: Could not open analysis tree file (" << name << ")" << endl;
        exit(-1);
      }
      sortData(inputfile, calfile, basisHPCoarse, basisHPFine, coarseBasisInd3DMap, basisInd3DMap, numSegHitsHist, pos3DMap, posXYMap, pos3DMap1Hit, pos3DMap2Hit, coarseChisqHist, fineChisqHist);
      inputfile->Close();
    }
    fclose(listFile);
  }else{
    TFile * inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      cout << "ERROR: Could not open analysis tree file (" << infile << ")" << endl;
      exit(-1);
    }
    sortData(inputfile, calfile, basisHPCoarse, basisHPFine, coarseBasisInd3DMap, basisInd3DMap, numSegHitsHist, pos3DMap, posXYMap, pos3DMap1Hit, pos3DMap2Hit, coarseChisqHist, fineChisqHist);
    inputfile->Close();
  }




  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  inp->Close();
  inp2->Close();
}

int main(int argc, char ** argv) {

  const char *afile, *basisfileCoarse, *basisfileFine, *outfile, *calfile;
  char *ext;

  if (argc < 2) {
    cout << endl << "This sortcode sorts interaction positions using a waveform basis (generated using the GammaTrackingMakeBasis code) using the grid search method." << endl << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasisGridSearch analysis_tree_file basis_file_coarse basis_file_fine cal_file output_file" << endl << endl;
    cout << "The analysis tree (containing the calibration data used to make the basis) is a required argument.  Omitting other arguments will cause the sortcode to fall back to default values." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    basisfileCoarse = "trackingWaveformBasisCoarse.root";
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  } else if (argc == 3) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = "trackingWaveformBasisFine.root";
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  } else if (argc == 4) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingBasisSortTestGridSearch.root";
  }else if (argc == 5) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = argv[4];
	  outfile = "trackingBasisSortTestGridSearch.root";
  } else if (argc == 6) {
	  afile = argv[1];
    basisfileCoarse = argv[2];
    basisfileFine = argv[3];
	  calfile = argv[4];
	  outfile = argv[5];
  } else if (argc > 6) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingSortFromBasis analysis_tree_file basis_file_coarse basis_file_fine cal_file output_file" << endl;
	  return 0;
  }

  cout << "Starting sortcode..." << endl;

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0) {
	  grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  cout << "Input file: " << afile << endl << "Coarse Waveform basis file: " << basisfileCoarse << endl << "Fine Waveform basis file: " << basisfileFine << endl;
  cout << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, false);
  }

  return 0;
}
