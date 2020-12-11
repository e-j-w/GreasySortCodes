//#include "TCanvas.h"
//#include "TApplication.h"
#include "common.h" //define all global variables here!


TApplication *theApp;

TH1D *coarseBasis[(Int_t)(VOXEL_BINS_R*COARSE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR*VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR)];
TH1D *fineBasis[(Int_t)(VOXEL_BINS_R*FINE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR)];
TH1I *numSegHitsHist, *coarseBasisIndMap, *basisIndMap;
TH1D *coarseChisqHist, *fineChisqHist;
TH2D *posXYMap, *posXYMap1Hit, *posXYMap2Hit, *coarseChisqEHist;
TH3D *pos3DMap, *pos3DMap1Hit, *pos3DMap2Hit;
TH3I *coarseBasisInd3DMap, *basisInd3DMap;
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
void fillAtHitIndBasis(const Int_t hitInd, const Int_t basisBinsR, const Int_t basisBinsAngle, const Int_t basisBinsZ, TH1I *basisIndMap, TH3I *basisInd3DMap){
  Int_t rBinFine, angleBinFine, zBinFine;
  if(hitInd >= 0){
    //a best-fit fine basis waveform was found, figure out the corresponding bin
    basisIndMap->Fill(hitInd);
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

void sortData(TFile *inputfile, const char *calfile, TH1I *basisHPCoarse, TH1I *basisHPFine){
  
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
      if((coreCharge <= BASIS_MIN_ENERGY)||(coreCharge <= 0)) continue; //bad energy
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
            //cout << "Entry " << jentry << ", multiple hits in one segment." << endl;
            goodWaveforms = false;
            break;
          }else{
            segsInData|=(one<<(tigress_hit->GetSegmentHit(i).GetSegment()-1));
          }
          if(tigress_hit->GetSegmentHit(i).GetWaveform()->size()!=SAMPLES){
            //cout << "Entry " << jentry << ", mismatched waveform sizes." << endl;
            goodWaveforms = false;
            break;
          }
          if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.5*coreCharge){
            evtSegHP|=(one<<(tigress_hit->GetSegmentHit(i).GetSegment()-1)); //GetSegment() is 1-indexed
            numSegHits++;
          }
        }
        if(evtSegHP==0){
          goodWaveforms = false;
        }
        //cout << "HP: " << evtSegHP << ", num seg hits: " << numSegHits << endl;
        //cout << "good: " << goodWaveforms << ", num seg hits: " << numSegHits << endl;
        if(goodWaveforms){

          Double_t minChisqCoarse = BIG_NUMBER;
          
          Int_t minIndCoarse[NSEG];
          for(int i=0;i<NSEG;i++){
            minIndCoarse[i] = -1;
          }

          if(numSegHits<=NSEG){  
            //cout << "entry " << jentry << ", 2 segment hits." << endl;
            Int_t hitInd[NSEG];
            for(int i=0;i<NSEG;i++){
              hitInd[i]=-1;
            }
            for(int i=0; i<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i++){
              if(coarseBasis[i]!=NULL){
                //scaled sum of multiple basis waveforms case
                //cout << "entry " << jentry << ", scaled sum case." << endl;
                Int_t hitSeg[NSEG];
                Double_t scaleFacHit[NSEG];
                
                hitSeg[0] = commonHitInHP(evtSegHP,(Int_t)basisHPCoarse->GetBinContent(i+1),NSEG);
                if(hitSeg[0]>=0){
                  bool hitsFound=true;
                  hitInd[0] = i;
                  for(int hitNum=1; hitNum<numSegHits; hitNum++){
                    for(int i2=i+1; i2<(coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ); i2++){
                      if(coarseBasis[i2]!=NULL){
                        hitSeg[hitNum] = commonHitInHP(evtSegHP,(Int_t)basisHPCoarse->GetBinContent(i2+1),NSEG);
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
                        Double_t wfrmSampleVal = fabs((segwf->at(k) - seg_waveform_baseline)/coreCharge);
                        Double_t basisSampleVal = 0.;
                        for(int hitNum=0; hitNum<numSegHits; hitNum++){
                          basisSampleVal += scaleFacHit[hitNum]*coarseBasis[hitInd[hitNum]]->GetBinContent(segNum*SAMPLES + k + 1);
                        }
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

                    //check if chisq is at minimum
                    if((chisq>0)&&(chisq<minChisqCoarse)){
                      minChisqCoarse = chisq;
                      for(int hitNum=0; hitNum<numSegHits; hitNum++){
                        minIndCoarse[hitNum] = hitInd[hitNum];
                      }
                    }

                  }

                }
              }
            }
            //cout << "Entry " << jentry << " chisq val: " << minChisqCoarse << ", min indices: " << minIndCoarse[0] << ", " << minIndCoarse[1] << endl;
          }else{
            cout << "WARNING: Entry " << jentry << " contains an unhandled number of hit segments (" << numSegHits << ")." << endl;
            continue;
          }
          //cout << "Entry " << jentry << ", num seg hits: " << numSegHits << ", coarse basis min chisq: " << minChisqCoarse << endl;

          //check that the coarse basis sort was successful
          bool hitsFound=true;
          for(int hitNum=0; hitNum<numSegHits; hitNum++){
            if(minIndCoarse[hitNum] < 0){
              //at least one hit didn't have a basis index found
              hitsFound = false;
              break;
            }
          }
          if(!hitsFound){
            continue;
          }
          
          coarseChisqHist->Fill(minChisqCoarse);
          coarseChisqEHist->Fill(minChisqCoarse,coreCharge);

          /*if(minChisqCoarse > 1){
            TCanvas *c1 = new TCanvas("c1","Histogram",200,10,1200,1000);
            TH1D *ep = new TH1D("event_pulse","event_pulse",SAMPLES*(NSEG),0,SAMPLES*(NSEG));
            for(int i = 0; i < NSEG; i++){
              Double_t seg_waveform_baseline = 0.;
              for(int k = 0; k < BASELINE_SAMPLES; k++){
                seg_waveform_baseline += tigress_hit->GetSegmentHit(i).GetWaveform()->at(k);
              }
              seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;
              for(int k = 0; k < SAMPLES; k++){
                ep->SetBinContent((tigress_hit->GetSegmentHit(i).GetSegment()-1)*SAMPLES + k,(tigress_hit->GetSegmentHit(i).GetWaveform()->at(k) - seg_waveform_baseline)/coreCharge);
              }
            }
            cout << "entry " << jentry << ", hit " << hitInd << ", event HP: " << evtSegHP << ", chisq: " << minChisqCoarse << ", core charge: " << coreCharge << endl;
            ep->Draw();
            theApp->Run(kTRUE);
          }*/

          //now do everything with the fine basis
          if(minChisqCoarse<MAX_BASIS_SORT_CHISQ){

            Int_t minIndFine[NSEG];
            for(int i=0;i<NSEG;i++){
              minIndFine[i] = -1;
            }
            Int_t hitInd[NSEG];
            for(int i=0;i<NSEG;i++){
              hitInd[i]=-1;
            }
            Int_t rBinCoarse[NSEG], angleBinCoarse[NSEG], zBinCoarse[NSEG];
            Int_t rBinFine[NSEG], angleBinFine[NSEG], zBinFine[NSEG], basisInd[NSEG];
            
            Double_t minChisqFine = BIG_NUMBER;

            for(int hitNum=0; hitNum<numSegHits; hitNum++){
              fillAtHitIndBasis(minIndCoarse[hitNum],coarseBasisBinsR,coarseBasisBinsAngle,coarseBasisBinsZ,coarseBasisIndMap,coarseBasisInd3DMap);

              //a best-fit coarse basis waveform was found, figure out the corresponding bin
              rBinCoarse[hitNum] = (Int_t)(minIndCoarse[hitNum]/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
              angleBinCoarse[hitNum] = (Int_t)((minIndCoarse[hitNum] - rBinCoarse[hitNum]*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
              zBinCoarse[hitNum] = (minIndCoarse[hitNum] - rBinCoarse[hitNum]*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBinCoarse[hitNum]*coarseBasisBinsZ);
              //rVal = (rBin+0.5)*MAX_VAL_R/(1.0*coarseBasisBinsR);
              //angleVal = (angleBin+0.5)*360./(1.0*coarseBasisBinsAngle);
              //zVal = (zBin+0.5)*MAX_VAL_Z/(1.0*coarseBasisBinsZ);
              /*if(rBinCoarse==0){
                cout << "coarse basis bins: [ " << rBinCoarse << " " << angleBinCoarse << " " << zBinCoarse << " ]" << endl;
                cout << "angle bins coarse: " << 4*getNumAngleBins(rBinCoarse,COARSE_BASIS_BINFACTOR) << endl;
              }*/

            }

            for(int i=0; i<(basisBinRatio*basisBinRatio*basisBinRatio); i++){

              rBinFine[0] = (rBinCoarse[0]*basisBinRatio) + (Int_t)(i/(basisBinRatio*basisBinRatio*1.0));
              angleBinFine[0] = (angleBinCoarse[0]*basisBinRatio) + (((Int_t)(i/(basisBinRatio*1.0)) % basisBinRatio) );
              zBinFine[0] = (zBinCoarse[0]*basisBinRatio) + (i % basisBinRatio);
              basisInd[0] = rBinFine[0]*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine[0]*fineBasisBinsZ + zBinFine[0];
              //cout << "coarse basis bins [" << rBinCoarse << " " << angleBinCoarse << " " << zBinCoarse << "]" << endl;
              //cout << "fine basis bin: " << basisInd << endl;
              /*if(rBinCoarse==0){
                cout << "fine basis bins [" << rBinFine << " " << angleBinFine << " " << zBinFine << "], angle bins fine : " << 4*getNumAngleBins(rBinFine,FINE_BASIS_BINFACTOR,COARSE_BASIS_BINFACTOR)*FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR << ", fine basis index: " << basisInd << endl;
              }*/

              if(fineBasis[basisInd[0]]!=NULL){
                //cout << "   fine basis index exists" << endl;
              
                //signal resprented by either a single basis waveform with both segments hit (typically near segment boundaries)
                //or by multiple basis waveforms with hits in different segments
                //in the latter case, should scale each basis waveform based on the corresponding segment energy and then 
                //add together before performing chisq test
              
                //scaled sum of multiple basis waveforms case
                Int_t hitSeg[NSEG];
                Double_t scaleFacHit[NSEG];

                hitSeg[0] = commonHitInHP(evtSegHP,(Int_t)basisHPFine->GetBinContent(basisInd[0]+1),NSEG);
                if(hitSeg[0]>=0){
                  bool hitsFound=true;
                  hitInd[0] = i;
                  for(int hitNum=1; hitNum<numSegHits; hitNum++){

                    rBinCoarse[hitNum] = (Int_t)(minIndCoarse[hitNum]/(coarseBasisBinsAngle*coarseBasisBinsZ*1.0));
                    angleBinCoarse[hitNum] = (Int_t)((minIndCoarse[hitNum] - rBinCoarse[hitNum]*(coarseBasisBinsAngle*coarseBasisBinsZ*1.0))/(coarseBasisBinsZ*1.0));
                    zBinCoarse[hitNum] = (minIndCoarse[hitNum] - rBinCoarse[hitNum]*(coarseBasisBinsAngle*coarseBasisBinsZ) - angleBinCoarse[hitNum]*coarseBasisBinsZ);

                    for(int i2=i+1; i2<(basisBinRatio*basisBinRatio*basisBinRatio); i2++){
                      
                      rBinFine[hitNum] = (rBinCoarse[hitNum]*basisBinRatio) + (Int_t)(i2/(basisBinRatio*basisBinRatio*1.0));
                      angleBinFine[hitNum] = (angleBinCoarse[hitNum]*basisBinRatio) + (((Int_t)(i2/(basisBinRatio*1.0)) % basisBinRatio) );
                      zBinFine[hitNum] = (zBinCoarse[hitNum]*basisBinRatio) + (i2 % basisBinRatio);
                      Int_t basisIndTmp = rBinFine[hitNum]*fineBasisBinsAngle*fineBasisBinsZ + angleBinFine[hitNum]*fineBasisBinsZ + zBinFine[hitNum];

                      if(fineBasis[basisIndTmp]!=NULL){

                        hitSeg[hitNum] = commonHitInHP(evtSegHP,(Int_t)basisHPFine->GetBinContent(basisIndTmp+1),NSEG);
                        if(hitSeg[hitNum]>=0){
                          bool goodSeg=true;
                          for(int k=0;k<hitNum;k++){
                            if(hitSeg[hitNum]==hitSeg[k]){
                              goodSeg=false;
                            }
                          }
                          if(goodSeg){
                            hitInd[hitNum] = i2;
                            basisInd[hitNum] = basisIndTmp;
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
                    }

                    //compute chisq
                    Double_t chisq = 0.;
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
                        Double_t wfrmSampleVal = fabs((segwf->at(k) - seg_waveform_baseline)/coreCharge);
                        Double_t basisSampleVal = 0.;
                        for(int hitNum=0; hitNum<numSegHits; hitNum++){
                          //cout << "here7, hitNum: " << hitNum << ", ind: " << basisInd[hitNum] << ", scalefac: " << scaleFacHit[hitNum] << endl;
                          basisSampleVal += scaleFacHit[hitNum]*fineBasis[basisInd[hitNum]]->GetBinContent(segNum*SAMPLES + k + 1);
                        }
                        chisq += pow(wfrmSampleVal - basisSampleVal,2)*hitWeight;
                      }
                    }

                    //check if chisq is at minimum
                    if((chisq>0)&&(chisq<minChisqFine)){
                      minChisqFine = chisq;
                      for(int hitNum=0; hitNum<numSegHits; hitNum++){
                        minIndFine[hitNum] = basisInd[hitNum];
                      }
                    }
                  }
                }
              }
            }
            /*if(numSegHits==2){
              cout << "2 hit event, fine basis indices: [" << minIndFine[0] << " " << minIndFine[1] << "]" << ", chisq: " << minChisqFine << endl;
            }*/

            fineChisqHist->Fill(minChisqFine);
            //fill histograms
            if(minChisqFine<MAX_BASIS_SORT_CHISQ){
              for(int i=0;i<2;i++){
                if(minIndFine[i]>=0){
                  numSegHitsHist->Fill(numSegHits);
                  fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap);
                  fillAtHitIndBasis(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,basisIndMap,basisInd3DMap);
                  fillXYAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,posXYMap);
                  if(numSegHits==1){
                    fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap1Hit);
                    fillXYAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,posXYMap1Hit);
                  }else if(numSegHits==2){
                    fillAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,randGen,pos3DMap2Hit);
                    fillXYAtHitInd(minIndFine[i],fineBasisBinsR,fineBasisBinsAngle,fineBasisBinsZ,posXYMap2Hit);
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

  numSegHitsHist = new TH1I("numSegHits","numSegHits",10,0,10);
  list->Add(numSegHitsHist);
  pos3DMap = new TH3D("pos3DMap","pos3DMap",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap);
  pos3DMap1Hit = new TH3D("pos3DMap1Hit","pos3DMap1Hit",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap1Hit);
  pos3DMap2Hit = new TH3D("pos3DMap2Hit","pos3DMap2Hit",40,-40,40,40,-40,40,40,-10,100);
  list->Add(pos3DMap2Hit);
  posXYMap = new TH2D("posXYMap","posXYMap",40,-40,40,40,-40,40);
  list->Add(posXYMap);
  posXYMap1Hit = new TH2D("posXYMap1Hit","posXYMap1Hit",40,-40,40,40,-40,40);
  list->Add(posXYMap1Hit);
  posXYMap2Hit = new TH2D("posXYMap2Hit","posXYMap2Hit",40,-40,40,40,-40,40);
  list->Add(posXYMap2Hit);
  coarseBasisIndMap = new TH1I("coarseBasisIndMap","coarseBasisIndMap",coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ,0,coarseBasisBinsR*coarseBasisBinsAngle*coarseBasisBinsZ);
  list->Add(coarseBasisIndMap);
  basisIndMap = new TH1I("basisIndMap","basisIndMap",fineBasisBinsR*fineBasisBinsAngle*fineBasisBinsZ,0,fineBasisBinsR*fineBasisBinsAngle*fineBasisBinsZ);
  list->Add(basisIndMap);
  coarseBasisInd3DMap = new TH3I("coarseBasisInd3DMap","coarseBasisInd3DMap",coarseBasisBinsR,0,coarseBasisBinsR,coarseBasisBinsAngle,0,coarseBasisBinsAngle,coarseBasisBinsZ,0,coarseBasisBinsZ);
  list->Add(coarseBasisInd3DMap);
  basisInd3DMap = new TH3I("basisInd3DMap","basisInd3DMap",fineBasisBinsR,0,fineBasisBinsR,fineBasisBinsAngle,0,fineBasisBinsAngle,fineBasisBinsZ,0,fineBasisBinsZ);
  list->Add(basisInd3DMap);
  coarseChisqHist = new TH1D("coarseChisqHist","coarseChisqHist",1000,0,200);
  list->Add(coarseChisqHist);
  fineChisqHist = new TH1D("fineChisqHist","fineChisqHist",1000,0,200);
  list->Add(fineChisqHist);
  coarseChisqEHist = new TH2D("coarseChisqEHist","coarseChisqEHist",1000,0,200,100,0,4000);
  list->Add(coarseChisqEHist);

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
      sortData(inputfile, calfile, basisHPCoarse, basisHPFine);
      inputfile->Close();
    }
    fclose(listFile);
  }else{
    TFile * inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      cout << "ERROR: Could not open analysis tree file (" << infile << ")" << endl;
      exit(-1);
    }
    sortData(inputfile, calfile, basisHPCoarse, basisHPFine);
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
    cout << "The analysis tree (containing the experimental data used to be sorted) is a required argument.  Omitting other arguments will cause the sortcode to fall ";
    cout << "back to default values.  The analysis tree can be a single ROOT file, or a list of ROOT files (using file extension .list) can be specified instead." << endl << endl;
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

  //theApp=new TApplication("App", &argc, argv);
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    sort_from_basis(afile, basisfileCoarse, basisfileFine, calfile, outfile, false);
  }

  return 0;
}
