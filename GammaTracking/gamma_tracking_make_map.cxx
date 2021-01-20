#include "GammaTrackingTIGRESS.h" //define all global variables here!

//forward declarations
Int_t getNumAngleBins(Int_t rInd, Double_t rScaleFac, Double_t scaleFac);
Int_t getNumAngleBins(Int_t rInd, Double_t scaleFac);
double calc_ordering(TTigressHit *, const Int_t, const Int_t);

char hname[64];
TH1D *rDistHist[NSEG], *angleDistHist[NSEG*VOXEL_BINS_R], *zDistHist[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];
TH1 *rDistHistC[NSEG], *angleDistHistC[NSEG*VOXEL_BINS_R], *zDistHistC[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];

void sortData(TFile *inputfile, const char *calfile, TH3D *rhophizetaHist[NSEG]){
  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  cout << AnalysisTree->GetNtrees() << " tree files, details:" << endl;
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  cout << "Reading calibration file: " << calfile << endl;
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TTigress *tigress = 0;
  TTigressHit *tigress_hit;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    exit(-1);
  }

  Int_t hit_counter = 0;
  Int_t map_hit_counter = 0;
  Int_t overflow_rho_counter = 0;
  Int_t overflow_phi_counter = 0;
  Int_t overflow_zeta_counter = 0;

  Int_t offset = 0;
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue;  //exclude pileup
      Double_t coreCharge = tigress_hit->GetCharge();
      if((coreCharge <= 0)||(coreCharge > BASIS_MAX_ENERGY)) continue; //bad energy
      hit_counter++; 
      //cout << "Number of segments: " << tigress_hit->GetSegmentMultiplicity() << endl;
      if(tigress_hit->GetSegmentMultiplicity() == NSEG){
        //all segments have waveforms      

        //check that the waveforms are the same size, and that there are no duplicate segments
        bool goodWaveforms = true;
        Int_t numSegHits = 0; //counter for the number of segments with a hit (ie. over the threshold energy)
        Int_t segsInData = 0;
        Int_t numSamples = -1;
        Double_t maxSegCharge = 0.;
        Int_t maxChargeSeg = -1;
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
          if((tigress_hit->GetSegmentHit(i).GetCharge() > BASIS_MAX_ENERGY)||(fabs(tigress_hit->GetSegmentHit(i).GetCharge()) > MAX_ENERGY_SINGLE_INTERACTION)){
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
        if(numSegHits != 1){
          goodWaveforms = false;
        }
        if(goodWaveforms){
            
          TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(maxChargeSeg);

          Int_t segNum = segment_hit.GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

          //calculate all ordering parameters (see ordering_parameter_calc.cxx)
          Double_t rho = calc_ordering(tigress_hit,maxChargeSeg,0);
          if(rho == BAD_RETURN){
            continue;
          }
          Double_t phi = calc_ordering(tigress_hit,maxChargeSeg,1);
          if(phi == BAD_RETURN){
            continue;
          }
          Double_t zeta = calc_ordering(tigress_hit,maxChargeSeg,2);
          if(zeta == BAD_RETURN){
            continue;
          }

          //cout << "seg: " << segNum << ", rho: " << rho << ", phi:" << phi << ", zeta: " << zeta << endl;

          if(fabs(rho) > RHO_MAX)
            overflow_rho_counter++;
          if(fabs(phi) > PHI_MAX)
            overflow_phi_counter++;
          if(fabs(zeta) > ZETA_MAX)
            overflow_zeta_counter++;
          if((fabs(rho) > RHO_MAX)||(fabs(phi) > PHI_MAX)||(fabs(zeta) > ZETA_MAX)){
            continue;
          }
          
          //cout << "bin: " << rhophizetaHist[segNum]->GetBin(rho,phi,zeta) << endl;
          rhophizetaHist[segNum]->Fill(rho,phi,zeta);
          map_hit_counter++;
          
        }
      }

      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << map_hit_counter << " of " << hit_counter << " hits retained (" << 100*map_hit_counter/hit_counter << " %)." << endl;
  cout << "Hits with ordering parameters out of map range for [rho, phi, zeta]: [" << overflow_rho_counter << " " << overflow_phi_counter << " " << overflow_zeta_counter << "]" << endl;
}


//function which generates a mapping between ordering parameters and real spatial coordinates
//and saves this mapping to disk
void generate_mapping(const char *infile, const char *simfile, const char *calfile, const char *outfile, bool inpList) {

  TList * list = new TList;
  
  //setup density functions for the real spatial coordinates r, angle, and z in cylindrical coordinates
  //these are read in from a ROOT tree generated by a GEANT4 simulation such as G4TIP
  for(int k = 0; k < NSEG; k++){
    sprintf(hname,"rDistSeg%i",k);
    rDistHist[k] = new TH1D(hname,Form("rDistSeg%i",k),160,0,MAX_VAL_R);
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      sprintf(hname,"angleDistSeg%ir%ito%i",k,rValMin,rValMax);
      angleDistHist[k*VOXEL_BINS_R + j] = new TH1D(hname,Form("angleDistSeg%ir%ito%i",k,rValMin,rValMax),120,0,MAX_VAL_ANGLE);
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i < numAngleBins; i++){
        Int_t angleValMin = (Int_t)(i*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        Int_t angleValMax = (Int_t)((i+1)*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        sprintf(hname,"zDistSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax);
        zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = new TH1D(hname,Form("zDistSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax),360,0,90);
      }
    }
  }

  TTree *simTree;
  TFile *inp = new TFile(simfile,"read");
  if((simTree = (TTree*)inp->Get("tree"))==NULL){
    cout << "ERROR: No spatial coordinate distribution info in the specified ROOT file!" << endl;
    exit(-1);
  }
  TBranch *rBranch, *phiBranch, *zBranch, *segIDBranch, *numHitsBranch;
  TLeaf *rLeaf, *phiLeaf, *zLeaf, *segIDLeaf, *numHitsLeaf;
  if((rBranch = simTree->GetBranch("TigressSegmentMaxECylSphr"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxECylr' doesn't correspond to a branch or leaf in the tree!" << endl;
    exit(-1);
  }else{
    rLeaf = (TLeaf*)rBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((phiBranch = simTree->GetBranch("TigressSegmentMaxECylSphphi"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxECylphi' doesn't correspond to a branch or leaf in the tree!" << endl;
    exit(-1);
  }else{
    phiLeaf = (TLeaf*)phiBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((zBranch = simTree->GetBranch("TigressSegmentMaxECylSphz"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxECylz' doesn't correspond to a branch or leaf in the tree!" << endl;
    exit(-1);
  }else{
    zLeaf = (TLeaf*)zBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((segIDBranch = simTree->GetBranch("TigressSegmentMaxEId"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxEId' doesn't correspond to a branch or leaf in the tree!" << endl;
    exit(-1);
  }else{
    segIDLeaf = (TLeaf*)segIDBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((numHitsBranch = simTree->GetBranch("TigressSegmentMaxENumHits"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxENumHits' doesn't correspond to a branch or leaf in the tree!" << endl;
    exit(-1);
  }else{
    numHitsLeaf = (TLeaf*)numHitsBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  
  Int_t sentries = simTree->GetEntries();
  Double_t rVal, angleVal, zVal;
  Int_t idVal, numHitsVal;
  for (int i=0;i<sentries;i++){
    simTree->GetEntry(i);
    for(int j=0; j<rLeaf->GetNdata(); j++) { //deal with multiple fold events
      if((rLeaf->GetNdata()==phiLeaf->GetNdata())&&((rLeaf->GetNdata()==zLeaf->GetNdata()))&&(rLeaf->GetNdata()==segIDLeaf->GetNdata())){
        rVal = rLeaf->GetValue(j);
        angleVal = phiLeaf->GetValue(j)*180./M_PI;
        zVal = zLeaf->GetValue(j);
        idVal = segIDLeaf->GetValue(j)-1; //convert to zero-indexed
        numHitsVal = numHitsLeaf->GetValue(j);
        //cout << "sim data seg ID: " << idVal << ", r: " << rVal << ", rInd: " << rInd << ", angle: " << angleVal << ", z: " << zVal << endl;
        if(numHitsVal==1){ //restrict to single interaction events
        //if(numHitsVal>0){
          if((idVal>=0)&&(idVal<NSEG)){
            if((rVal==rVal)&&(rVal>=0)&&(rVal<=MAX_VAL_R)){ 
              rDistHist[idVal]->Fill(rVal);
              Int_t rInd = (Int_t)(rVal*VOXEL_BINS_R/MAX_VAL_R);
              if((angleVal==angleVal)&&(angleVal>=0)&&(angleVal<=MAX_VAL_ANGLE)){
                angleDistHist[idVal*VOXEL_BINS_R + rInd]->Fill(angleVal);
                Int_t angleInd = (Int_t)(angleVal*getNumAngleBins(rInd,1.0)/MAX_VAL_ANGLE);
                zDistHist[idVal*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]->Fill(zVal);
              }
            }
          }else{
            cout << "Sim tree: invalid segment." << endl;
          }
        }
      }else{
        cout << "Sim tree: data size mismatch." << endl;
      }

    }
  }

  //normalize distributions and get cumulative versions
  
  for(int k = 0; k < NSEG; k++){
    rDistHist[k]->Scale(1.0/rDistHist[k]->Integral());
    rDistHistC[k] = rDistHist[k]->GetCumulative();
    list->Add(rDistHist[k]);
    for(int j = 0; j < VOXEL_BINS_R; j++){
      //printf("ind: %i, integral: %f\n",angleDistHist[k*VOXEL_BINS_R + j]->Integral())
      if(angleDistHist[k*VOXEL_BINS_R + j]->Integral()>0.){
        angleDistHist[k*VOXEL_BINS_R + j]->Scale(1.0/angleDistHist[k*VOXEL_BINS_R + j]->Integral());
        angleDistHistC[k*VOXEL_BINS_R + j] = angleDistHist[k*VOXEL_BINS_R + j]->GetCumulative();
        list->Add(angleDistHist[k*VOXEL_BINS_R + j]);
      }
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i < numAngleBins; i++){
        if(zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->Integral()>0.){
          zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->Scale(1.0/zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->Integral());
          zDistHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetCumulative();
          list->Add(zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]);
        }
        
      }
    }
  }

  cout << "Spatial parameter distributions read in." << endl;

  //setup a histogram for the ordering parameters
  //ROOT histograms are used to store this data since ROOT provides useful
  //methods such as GetCumulative() and GetQuantiles() which will be used later
  //TH3D is used to preverse correlated values of rho, phi, zeta, 
  //will handle the individual mapping of the parameters later
  TH3D *rhophizetaHist[NSEG];
  for(int k = 0; k < NSEG; k++){
    sprintf(hname,"rhophizetaSeg%i",k);
    rhophizetaHist[k] = new TH3D(hname,Form("rhophizetaSeg%i",k),N_BINS_ORDERING,-1.*RHO_MAX,1.*RHO_MAX,N_BINS_ORDERING,-1.*PHI_MAX,1.*PHI_MAX,N_BINS_ORDERING,-1.*ZETA_MAX,1.*ZETA_MAX);
    //list->Add(rhophizetaHist[k]);
  }


  if(inpList){
    FILE *listFile;
    char name[256];
    if((listFile=fopen(infile,"r"))==NULL){
      cout << "ERROR: Could not open analysis tree list file!" << endl;
      exit(-1);
    }
    while(fscanf(listFile,"%s",name)!=EOF){
      TFile *inputfile = new TFile(name, "READ");
      if (!inputfile->IsOpen()) {
        cout << "ERROR: Could not open analysis tree file (" << name << ")" << endl;
        exit(-1);
      }
      sortData(inputfile,calfile,rhophizetaHist);
      inputfile->Close();
    }
    fclose(listFile);
  }else{
    TFile *inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      cout << "ERROR: Could not open analysis tree file (" << infile << ")" << endl;
      exit(-1);
    }
    sortData(inputfile,calfile,rhophizetaHist);
    inputfile->Close();
  }
  

  //generate normalized cumulative distributions of all ordering parameters,
  //and then ordering parameter to spatial coordinate maps
  Double_t xVal[1], qVal[1];
  Int_t nQuantiles;

  //first we map rho on a per segment basis
  cout << "Mapping ordering parameter rho to radius spatial parameter..." << endl;
  //generate cumulative distributions
  TH1 *rhoHist[NSEG], *rhoHistC[NSEG];
  for(int k = 0; k < NSEG; k++){
    rhoHist[k] = new TH1D();
    rhoHist[k] = rhophizetaHist[k]->ProjectionX();
    if(rhoHist[k]->GetEntries()>0){
      rhoHist[k]->Scale(1.0/rhoHist[k]->Integral());
      rhoHistC[k] = rhoHist[k]->GetCumulative();
      rhoHistC[k]->SetNameTitle(Form("rhoHistCumulativeSeg%i",k),Form("rhoHistCumulativeSeg%i",k));
      list->Add(rhoHistC[k]);
    }
  }
  //generate maps
  TH1 *rMap[NSEG];
  for(int k = 0; k < NSEG; k++){
    if(rhoHist[k]->GetEntries()>0){
      sprintf(hname,"rMapSeg%i",k);
      rMap[k] = new TH1D(hname,Form("rMapSeg%i",k),N_BINS_ORDERING,-1.*RHO_MAX,1.*RHO_MAX);
      for(int l=0;l<N_BINS_ORDERING;l++){
        xVal[0] = rhoHistC[k]->GetBinContent(l+1);
        nQuantiles = rDistHist[k]->GetQuantiles(1,qVal,xVal);
        if(nQuantiles==1){
          rMap[k]->SetBinContent(l+1,qVal[0]);
        }
      }
      list->Add(rMap[k]);
    }
  }
  //store rho bin numbers which correspond to discrete values of r
  Int_t rhoBinVal[NSEG*(VOXEL_BINS_R + 1)];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j <= VOXEL_BINS_R; j++){
      Double_t rBinVal = j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)); //keep r bins equal in size
      rhoBinVal[k*(VOXEL_BINS_R+1) + j] = rMap[k]->FindFirstBinAbove(rBinVal-0.01); //-0.01 to handle bin values at the end of the mapping range
      //get correct range for last bin (which may extend radially beyond the physical detector edge)
      if((j > 0)&&(rhoBinVal[k*(VOXEL_BINS_R+1) + j] == -1)){
        rhoBinVal[k*(VOXEL_BINS_R+1) + j] = rMap[k]->FindLastBinAbove((j-1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)) - 0.01);
      }
      
    }
    /*for(int j = 0; j <= VOXEL_BINS_R; j++){
      cout << "seg: " << k << ", rhoBin: " << j <<  ", rhoBinVal: " << rhoBinVal[k*(VOXEL_BINS_R+1) + j] << endl;
    }*/
  }

  //now map phi based on mapped radius (binned) and segment number
  cout << "Mapping ordering parameter phi to angle spatial parameter..." << endl;
  //generate cumulative distributions
  TH1 *phiHist[NSEG*VOXEL_BINS_R], *phiHistC[NSEG*VOXEL_BINS_R];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j < VOXEL_BINS_R; j++){
      //get the projection in phi (y) gated on a range in rho (x)
      phiHist[k*VOXEL_BINS_R + j] = new TH1D();
      phiHistC[k*VOXEL_BINS_R + j] = new TH1D();
      Int_t lowerRhoBound = rhoBinVal[k*(VOXEL_BINS_R+1) + j];
      Int_t upperRhoBound = rhoBinVal[k*(VOXEL_BINS_R+1) + j + 1];
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      //cout << "seg: " << k << ", bin: " << j << ", proj bounds: " << lowerRhoBound << " " << upperRhoBound << endl;
      if((lowerRhoBound > 0)&&(lowerRhoBound < upperRhoBound)){
        phiHist[k*VOXEL_BINS_R + j] = rhophizetaHist[k]->ProjectionY("",lowerRhoBound,upperRhoBound);
        if(phiHist[k*VOXEL_BINS_R + j]->GetEntries()>0){
          phiHist[k*VOXEL_BINS_R + j]->Scale(1.0/phiHist[k*VOXEL_BINS_R + j]->Integral());
          phiHistC[k*VOXEL_BINS_R + j] = phiHist[k*VOXEL_BINS_R + j]->GetCumulative();
          phiHistC[k*VOXEL_BINS_R + j]->SetNameTitle(Form("phiHistCumulativeSeg%ir%ito%i",k,rValMin,rValMax),Form("phiHistCumulativeSeg%ir%ito%i",k,rValMin,rValMax));
          list->Add(phiHistC[k*VOXEL_BINS_R + j]);
        }
      }
    }
  }
  //generate maps
  TH1 *angleMap[NSEG*VOXEL_BINS_R];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      sprintf(hname,"angleMapSeg%ir%ito%i",k,rValMin,rValMax);
      angleMap[k*VOXEL_BINS_R + j] = new TH1D(hname,Form("angleMapSeg%ir%ito%i",k,rValMin,rValMax),N_BINS_ORDERING,-1.0*PHI_MAX,PHI_MAX);
      if((phiHist[k*VOXEL_BINS_R + j]->GetEntries()>0)&&(angleDistHist[k*VOXEL_BINS_R + j]->Integral()>0.)){
        for(int l=0;l<N_BINS_ORDERING;l++){
          xVal[0] = phiHistC[k*VOXEL_BINS_R + j]->GetBinContent(l+1);
          nQuantiles = angleDistHist[k*VOXEL_BINS_R + j]->GetQuantiles(1,qVal,xVal);
          if(nQuantiles>=1){
            //printf("val=%f ",qVal[0]); //=0 when there is no GEANT4 data in this bin
            angleMap[k*VOXEL_BINS_R + j]->SetBinContent(l+1,qVal[0]);
          }
        }
        list->Add(angleMap[k*VOXEL_BINS_R + j]);
      }
    }
  }
  //store bin numbers which correspond to discrete values of the angle
  Int_t phiBinVal[NSEG*(VOXEL_BINS_ANGLE_MAX + 1)*(VOXEL_BINS_R)];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i <= numAngleBins; i++){
        //printf("k j i %i %i %i\n",k,j,i);
        //printf("bins %i %i\n",k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i, k*VOXEL_BINS_R + j);
        Double_t angleBinVal = i*(MAX_VAL_ANGLE/(numAngleBins*1.0)); //size of phi bins depends on r
        phiBinVal[k*(VOXEL_BINS_ANGLE_MAX+1)*(VOXEL_BINS_R) + j*(VOXEL_BINS_ANGLE_MAX+1) + i] = angleMap[k*VOXEL_BINS_R + j]->FindFirstBinAbove(angleBinVal-0.01); //-0.01 to handle bin values at the end of the mapping range
        //printf("k j i %i %i %i, binVal: %i\n",k,j,i,phiBinVal[k*(VOXEL_BINS_ANGLE_MAX+1)*(VOXEL_BINS_R) + j*(VOXEL_BINS_ANGLE_MAX+1) + i]);
      }
    }
  }

  //now map zeta based on mapped radius and angle (both binned) and segment number
  cout << "Mapping ordering parameter zeta to depth spatial parameter..." << endl;
  //generate cumulative distributions
  TH1 *zetaHist[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)], *zetaHistC[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t lowerRhoBound = rhoBinVal[k*(VOXEL_BINS_R+1) + j];
      Int_t upperRhoBound = rhoBinVal[k*(VOXEL_BINS_R+1) + j + 1];
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i < numAngleBins; i++){
        //get the projection in zeta (z) gated on a range in rho (x) and phi (y)
        zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = new TH1D();
        zetaHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = new TH1D();
        Int_t lowerPhiBound = phiBinVal[k*(VOXEL_BINS_ANGLE_MAX+1)*(VOXEL_BINS_R) + j*(VOXEL_BINS_ANGLE_MAX+1) + i];
        Int_t upperPhiBound = phiBinVal[k*(VOXEL_BINS_ANGLE_MAX+1)*(VOXEL_BINS_R) + j*(VOXEL_BINS_ANGLE_MAX+1) + i + 1];
        Int_t angleValMin = (Int_t)(i*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        Int_t angleValMax = (Int_t)((i+1)*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        //cout << "bin: " << k << " " << j << " " << i << ", proj bounds x: " << lowerBoundx << " " << upperBoundx << endl;
        //cout << "proj bounds y: " << lowerBoundy << " " << upperBoundy << endl;
        if((lowerRhoBound > 0)&&(lowerRhoBound < upperRhoBound)){
          if((lowerPhiBound > 0)&&(lowerPhiBound < upperPhiBound)){
            zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = rhophizetaHist[k]->ProjectionZ("",lowerRhoBound,upperRhoBound,lowerPhiBound,upperPhiBound);
            if(zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetEntries()>0){
              zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->Scale(1.0/zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX]->Integral());
              zetaHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetCumulative();
              zetaHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->SetNameTitle(Form("zetaHistCumulativeSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax),Form("zetaHistCumulativeSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax));
              list->Add(zetaHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]);
            }
          }
        }
      }
    }
  }
  //generate maps
  TH1 *zMap[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];
  for(int k = 0; k < NSEG; k++){
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i < numAngleBins; i++){
        //cout << "bin: " << k << " " << j << " " << i << endl;
        Int_t angleValMin = (Int_t)(i*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        Int_t angleValMax = (Int_t)((i+1)*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        if((zetaHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetEntries()>0)&&(zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->Integral()>0.)){
          sprintf(hname,"zMapSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax);
          zMap[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = new TH1D(hname,Form("zMapSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax),N_BINS_ORDERING,-1.*ZETA_MAX,1.*ZETA_MAX);
          for(int l=0;l<N_BINS_ORDERING;l++){
            xVal[0] = zetaHistC[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetBinContent(l+1);
            nQuantiles = zDistHist[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->GetQuantiles(1,qVal,xVal);
            if(nQuantiles>=1){
              zMap[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]->SetBinContent(l+1,qVal[0]);
            }
          }
          list->Add(zMap[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i]);
        }
      }
    }
  }

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, mapping complete!" << endl;
  myfile->Close();
  inp->Close();
}

int main(int argc, char ** argv) {

  const char *afile, *simfile, *outfile, *calfile;
  char *ext;

  if (argc < 3) {
    cout << endl << "This sortcode makes a map for the gamma tracking direct method which transforms ordering parameters to real spatial coordinates.  ";
    cout << "These maps are constructed from real calibration data (which should be provided in the analysis tree, the first argument), and simulated distributions of the ";
    cout << "spatial coordinates (which should be provided in a ROOT tree, the second argument).  GEANT4 simulations exist which can produce the simulated distribution data ";
    cout << "(for example, G4TIP (https://github.com/e-j-w/G4TIP/))." << endl << endl;
    cout << "Arguments: ./GammaTrackingMakeMap analysis_tree_file sim_tree_file cal_file output_file" << endl << endl;
    cout << "The analysis tree and simulation tree are required arguments.  Omitting other arguments will cause the sortcode to fall back to default values. The analyisis ";
    cout << "tree can be a single ROOT file, or a list of ROOT files (using file extension .list) can be specified instead." << endl << endl;
    cout << "NOTE: this code requires a LOT of memory (around 10GB with N_BINS_ORDERING set to 512 in common.h, scaling as N_BINS_ORDERING^3) to run." << endl << endl;
	  return 0;
  } else if (argc == 3) {
	  afile = argv[1];
    simfile = argv[2];
	  calfile = "CalibrationFile.cal";
	  outfile = "trackingMap.root";
  } else if (argc == 4) {
	  afile = argv[1];
    simfile = argv[2];
	  calfile = argv[3];
	  outfile = "trackingMap.root";
  } else if (argc == 5) {
	  afile = argv[1];
    simfile = argv[2];
	  calfile = argv[3];
	  outfile = argv[4];
  } else if (argc > 5) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingMakeMap analysis_tree_file sim_tree_file cal_file output_file" << endl;
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

  cout << "Input file: " << afile << endl << "Simulation data file: " << simfile << endl << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    generate_mapping(afile, simfile, calfile, outfile, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    generate_mapping(afile, simfile, calfile, outfile, false);
  }
  

  return 0;
}
