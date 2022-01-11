#include "GammaTrackingTIGRESS.h" //define all global variables here!

//forward declarations
Int_t getNumAngleBins(Int_t rInd, Double_t rScaleFac, Double_t scaleFac);
Int_t getNumAngleBins(Int_t rInd, Double_t scaleFac);
double calc_ordering(TTigressHit *, const Int_t, const Int_t);

char hname[64];

void sortData(TFile *inputfile, const char *calfile, TH3D *rhophizetaHist[NPOS*NCORE*NSEG]){
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

  Long64_t hit_counter = 0;
  Long64_t map_hit_counter = 0;
  Long64_t overflow_rho_counter = 0;
  Long64_t overflow_phi_counter = 0;
  Long64_t overflow_zeta_counter = 0;

  for(Long64_t jentry = 0; jentry < tree->GetEntries(); jentry++){
    tree->GetEntry(jentry);
    for(int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++){
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
          if(numSamples != (Int_t)tigress_hit->GetSegmentHit(i).GetWaveform()->size()){
            //cout << "Entry " << jentry << ", mismatched waveform size (" << tigress_hit->GetSegmentHit(i).GetWaveform()->size() << ")." << endl;
            goodWaveforms = false;
            break;
          }
          if((tigress_hit->GetSegmentHit(i).GetCharge() > BASIS_MAX_ENERGY)||(fabs(tigress_hit->GetSegmentHit(i).GetCharge()) > MAX_ENERGY_SINGLE_INTERACTION)){
            //cout << "Entry " << jentry << ", charge out of range." << endl;
            goodWaveforms = false;
            break;
          }
          if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.3*coreCharge){
            if(tigress_hit->GetSegmentHit(i).GetCharge() > maxSegCharge){
              maxSegCharge = tigress_hit->GetSegmentHit(i).GetCharge();
              maxChargeSeg = i;
            }
            numSegHits++;
          }
        }
        if((goodWaveforms)&&(numSegHits != 1)){
          //cout << "Entry " << jentry << ", incorrect number of segments hit (" << numSegHits << ")." << endl;
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

          Int_t arrayPos = tigress_hit->GetArrayNumber();
          if((arrayPos>=0)&&(arrayPos<(NPOS*NCORE))){
            //cout << "bin: " << rhophizetaHist[arrayPos*NSEG + segNum]->GetBin(rho,phi,zeta) << endl;
            rhophizetaHist[arrayPos*NSEG + segNum]->Fill(rho,phi,zeta);
            map_hit_counter++;
          }
          
          
          
        }
      }

      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << map_hit_counter << " of " << hit_counter << " hits retained (" << 100*map_hit_counter/hit_counter << " %)." << endl;
  cout << "Hits with ordering parameters out of map range for [rho, phi, zeta]: [" << overflow_rho_counter << " " << overflow_phi_counter << " " << overflow_zeta_counter << "]" << endl;
}

void read_sim_file(const char *simfile, TH1D *zDistHist[NPOS*NCORE*NSEG], TH1D *rDistHist[NPOS*NCORE*NSEG*VOXEL_BINS_Z], TH1D *angleDistHist[NPOS*NCORE*NSEG*(VOXEL_BINS_Z)*(VOXEL_BINS_R)]){
  cout << "Reading simulation tree in file: " << simfile << endl;
  TTree *simTree;
  TFile *inp = new TFile(simfile,"read");
  if((simTree = (TTree*)inp->Get("tree"))==NULL){
    cout << "ERROR: No spatial coordinate distribution info in the specified ROOT file (" << simfile << ")!" << endl;
    exit(-1);
  }
  TBranch *rBranch, *phiBranch, *zBranch, *posIDBranch, *coreIDBranch, *segIDBranch, *numHitsBranch;
  TLeaf *rLeaf, *phiLeaf, *zLeaf, *posIDLeaf, *coreIDLeaf, *segIDLeaf, *numHitsLeaf;
  if((rBranch = simTree->GetBranch("TigressSegmentFullECylr"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullECylr' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    rLeaf = (TLeaf*)rBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((phiBranch = simTree->GetBranch("TigressSegmentFullECylphi"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullECylphi' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    phiLeaf = (TLeaf*)phiBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((zBranch = simTree->GetBranch("TigressSegmentFullECylz"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentMaxECylz' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    zLeaf = (TLeaf*)zBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((posIDBranch = simTree->GetBranch("TigressSegmentFullEPosId"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullEPosId' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    posIDLeaf = (TLeaf*)posIDBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((coreIDBranch = simTree->GetBranch("TigressSegmentFullECoreId"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullECoreId' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    coreIDLeaf = (TLeaf*)coreIDBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((segIDBranch = simTree->GetBranch("TigressSegmentFullESegId"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullESegId' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    segIDLeaf = (TLeaf*)segIDBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  if((numHitsBranch = simTree->GetBranch("TigressSegmentFullENumHits"))==NULL){
    cout << "ERROR: Sort data path 'TigressSegmentFullENumHits' doesn't correspond to a branch or leaf in the tree (" << simfile << ")!" << endl;
    exit(-1);
  }else{
    numHitsLeaf = (TLeaf*)numHitsBranch->GetListOfLeaves()->First(); //get the first leaf from the specified branch       
  }
  Long64_t sentries = simTree->GetEntries();
  Double_t rVal, angleVal, zVal;
  Int_t arrayPosVal, segIDVal, numHitsVal;
  for(Long64_t i=0;i<sentries;i++){
    simTree->GetEntry(i);
    for(Long64_t j=0; j<rLeaf->GetNdata(); j++) { //deal with multiple fold events
      if((rLeaf->GetNdata()==phiLeaf->GetNdata())&&((rLeaf->GetNdata()==zLeaf->GetNdata()))&&(rLeaf->GetNdata()==segIDLeaf->GetNdata())&&(rLeaf->GetNdata()==coreIDLeaf->GetNdata())&&(rLeaf->GetNdata()==posIDLeaf->GetNdata())&&(rLeaf->GetNdata()==numHitsLeaf->GetNdata())){
        rVal = rLeaf->GetValue(j);
        angleVal = phiLeaf->GetValue(j)*180./M_PI;
        if(angleVal >= MAX_VAL_ANGLE){
          angleVal = MAX_VAL_ANGLE - 0.01;
        }
        zVal = zLeaf->GetValue(j);
        if(zVal >= MAX_VAL_Z){
          zVal =  0.01;
        }
        //transform r to value along E-field lines
        //cout << "rVal start: " << rVal << ", z: " << zVal << endl;
        rVal = getREFieldFromR(rVal,zVal);
        //cout << "rVal EF: " << rVal << ", transformed back: " << getRFromREField(rVal,zVal) << endl;
        if(rVal >= MAX_VAL_R){
          //cout << "rVal (transformed): " << rVal << endl;
          rVal = MAX_VAL_R - 0.01;
        }else if(rVal < 0.){
          //cout << "rVal (transformed): " << rVal << endl;
          rVal = 0.;
        }
        arrayPosVal = (posIDLeaf->GetValue(j)-1)*NCORE + (coreIDLeaf->GetValue(j)-1);
        segIDVal = segIDLeaf->GetValue(j)-1; //convert to zero-indexed
        numHitsVal = numHitsLeaf->GetValue(j);
        //cout << "sim data seg ID: " << segIDVal << ", r: " << rVal << ", angle: " << angleVal << ", z: " << zVal << endl;
        //if(numHitsVal==1){ //restrict to single interaction events
        if(numHitsVal>0){
          if((arrayPosVal>=0)&&(arrayPosVal<(NPOS*NCORE))){
            if((segIDVal>=0)&&(segIDVal<NSEG)){
              if((zVal==zVal)&&(zVal>=0)&&(zVal<MAX_VAL_Z)){ 
                zDistHist[arrayPosVal*NSEG + segIDVal]->Fill(zVal);
                Int_t zInd = (Int_t)(zVal*VOXEL_BINS_Z/MAX_VAL_Z);
                if((rVal==rVal)&&(rVal>=0)&&(rVal<MAX_VAL_R)){
                  rDistHist[arrayPosVal*NSEG*VOXEL_BINS_Z + segIDVal*VOXEL_BINS_Z + zInd]->Fill(rVal);
                  Int_t rInd = (Int_t)(rVal*VOXEL_BINS_R/MAX_VAL_R);
                  if((angleVal==angleVal)&&(angleVal>=0)&&(angleVal<MAX_VAL_ANGLE)){
                    angleDistHist[arrayPosVal*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + segIDVal*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + zInd*VOXEL_BINS_R + rInd]->Fill(angleVal);
                  }
                }
              }
            }else{
              cout << "Sim tree (" << simfile << "): invalid segment." << endl;
            }
          }else{
            cout << "Sim tree (" << simfile << "): invalid array position." << endl;
          }
        }
      }else{
        cout << "Sim tree (" << simfile << "): data size mismatch." << endl;
      }

    }
  }
  inp->Close();
}



//function which generates a mapping between ordering parameters and real spatial coordinates
//and saves this mapping to disk
void generate_mapping(const char *infile, const char *simfile, const char *calfile, const char *outfile, bool inpList, bool simList) {

  //allocate very large arrays of object pointers
  TH1D **zDistHist = new TH1D*[NPOS*NCORE*NSEG];
  TH1D **rDistHist = new TH1D*[NPOS*NCORE*NSEG*VOXEL_BINS_Z];
  TH1D **angleDistHist = new TH1D*[NPOS*NCORE*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z)];
  TH1 **zDistHistC = new TH1*[NPOS*NCORE*NSEG];
  TH1 **rDistHistC = new TH1*[NPOS*NCORE*NSEG*VOXEL_BINS_Z];
  TH1 **angleDistHistC = new TH1*[NPOS*NCORE*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z)];

  TList *writeList = new TList;
  
  //setup density functions for the real spatial coordinates r, angle, and z in cylindrical coordinates
  //these are read in from a ROOT tree generated by a GEANT4 simulation such as G4TIP
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      sprintf(hname,"rDistPos%iSeg%i",l,k);
      zDistHist[l*(NSEG) + k] = new TH1D(hname,Form("zDistPos%iSeg%i",l,k),100,0,MAX_VAL_Z);
      for(int j = 0; j < VOXEL_BINS_Z; j++){
        Int_t zValMin = (Int_t)(j*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        Int_t zValMax = (Int_t)((j+1)*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        sprintf(hname,"rDistPos%iSeg%iz%ito%i",l,k,zValMin,zValMax);
        rDistHist[l*(NSEG)*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = new TH1D(hname,Form("rDistPos%iSeg%iz%ito%i",l,k,zValMin,zValMax),100,0,MAX_VAL_R);
        for(int i = 0; i < VOXEL_BINS_R; i++){
          Int_t rValMin = (Int_t)(i*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
          Int_t rValMax = (Int_t)((i+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
          sprintf(hname,"angleDistPos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax);
          angleDistHist[l*(NSEG)*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = new TH1D(hname,Form("angleDistPos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax),100,0,MAX_VAL_ANGLE);
        }
      }
    }
  }
  
  if(simList){
    FILE *listFile;
    char name[256];
    if((listFile=fopen(simfile,"r"))==NULL){
      cout << "ERROR: Could not open simulation tree list file!" << endl;
      exit(-1);
    }
    while(fscanf(listFile,"%s",name)!=EOF){
      read_sim_file(name,zDistHist,rDistHist,angleDistHist);
    }
    fclose(listFile);
  }else{
    read_sim_file(simfile,zDistHist,rDistHist,angleDistHist);
  }

  //normalize distributions and get cumulative versions
  cout << "Generating normalized and cumulative spatial parameter distributions..." << endl;
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      zDistHist[l*NSEG + k]->Scale(1.0/zDistHist[l*NSEG + k]->Integral());
      zDistHistC[l*NSEG + k] = zDistHist[l*NSEG + k]->GetCumulative();
      //writeList->Add(rDistHist[l*NSEG + k]);
      for(int j = 0; j < VOXEL_BINS_Z; j++){
        //printf("ind: %i, integral: %f\n",rDistHist[k*VOXEL_BINS_Z + j]->Integral())
        if(rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Integral()>0.){
          rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Scale(1.0/rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Integral());
          rDistHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetCumulative();
          //writeList->Add(rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]);
        }
        for(int i = 0; i < VOXEL_BINS_R; i++){
          if(angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Integral()>0.){
            angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Scale(1.0/angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Integral());
            angleDistHistC[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->GetCumulative();
            //writeList->Add(angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]);
          }
          
        }
      }
    }
  }
  
  cout << "Spatial parameter distributions read in." << endl;
  //setup histograms for the ordering parameters
  //ROOT histograms are used to store this data since ROOT provides useful
  //methods such as GetCumulative() and GetQuantiles() which will be used later
  //TH3D is used to preverse correlated values of rho, phi, zeta, 
  //will handle the individual mapping of the parameters later
  TH3D *rhophizetaHist[NPOS*NCORE*NSEG];
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      sprintf(hname,"rhophizetaPos%iSeg%i",l,k);
      rhophizetaHist[l*NSEG + k] = new TH3D(hname,Form("rhophizetaPos%iSeg%i",l,k),N_BINS_ORDERING,-1.*RHO_MAX,1.*RHO_MAX,N_BINS_ORDERING,-1.*PHI_MAX,1.*PHI_MAX,N_BINS_ORDERING,-1.*ZETA_MAX,1.*ZETA_MAX);
      //writeList->Add(rhophizetaHist[k]);
    }
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

  //first we map zeta on a per segment basis
  cout << "Mapping ordering parameter zeta to radius spatial parameter..." << endl;
  //generate cumulative distributions
  TH1 **zetaHist = new TH1*[NPOS*NCORE*NSEG];
  TH1 **zetaHistC = new TH1*[NPOS*NCORE*NSEG];
  TH1 **zMap = new TH1*[NPOS*NCORE*NSEG];
  Int_t *zetaBinVal = (Int_t*)malloc((NPOS*NCORE*NSEG*(VOXEL_BINS_Z + 1))*sizeof(Int_t));
  memset(zetaBinVal,0,sizeof(&zetaBinVal));
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      zetaHist[l*NSEG + k] = new TH1D();
      zetaHistC[l*NSEG + k] = new TH1D();
      zetaHist[l*NSEG + k] = rhophizetaHist[l*NSEG + k]->ProjectionZ("",1,N_BINS_ORDERING,1,N_BINS_ORDERING);
      if(zetaHist[l*NSEG + k]->GetEntries()>0){
        //generate cumulative distributions
        //cout << l << " entries: " << zetaHist[l*NSEG + k]->GetEntries() << ", integral: " << zetaHist[l*NSEG + k]->Integral() << endl;
        zetaHist[l*NSEG + k]->Scale(1.0/zetaHist[l*NSEG + k]->Integral());
        //cout << "entries2: " << zetaHist[l*NSEG + k]->GetEntries() << ", integral2: " << zetaHist[l*NSEG + k]->Integral() << endl;;
        zetaHistC[l*NSEG + k] = zetaHist[l*NSEG + k]->GetCumulative();
        //zetaHistC[l*NSEG + k]->SetNameTitle(Form("zetaHistCumulativePos%iSeg%i",l,k),Form("zetaHistCumulativePos%iSeg%i",l,k));
        //writeList->Add(zetaHistC[l*NSEG + k]);
        //generate maps
        sprintf(hname,"zMapPos%iSeg%i",l,k);
        zMap[l*NSEG + k] = new TH1D(hname,Form("zMapPos%iSeg%i",l,k),N_BINS_ORDERING,-1.*ZETA_MAX,1.*ZETA_MAX);
        for(int m=0;m<N_BINS_ORDERING;m++){
          xVal[0] = zetaHistC[l*NSEG + k]->GetBinContent(m+1);
          nQuantiles = zDistHist[l*NSEG + k]->GetQuantiles(1,qVal,xVal);
          if(nQuantiles==1){
            //cout << "m = " << m << ", nQuantiles = " << nQuantiles << ", val = " << qVal[0] << endl; //qVal[0]=0 when there is no GEANT4 data in this bin
            zMap[l*NSEG + k]->SetBinContent(m+1,qVal[0]);
          }
        }
        //cout << "here write zMap" << l*NSEG + k << endl;
        writeList->Add(zMap[l*NSEG + k]);
        //store zeta bin numbers which correspond to discrete values of z
        for(int j = 0; j <= VOXEL_BINS_Z; j++){
          Double_t zBinVal = j*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)); //keep z bins equal in size
          zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j] = zMap[l*NSEG + k]->FindFirstBinAbove(zBinVal-0.01); //-0.01 to handle bin values at the end of the mapping range
          //get correct range for last bin (which may extend radially beyond the physical detector edge)
          if((j > 0)&&(zetaBinVal[k*(VOXEL_BINS_Z+1) + j] == -1)){
            zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j] = zMap[l*NSEG + k]->FindLastBinAbove((j-1)*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)) - 0.01);
          }
        }
        /*for(int j = 0; j <= VOXEL_BINS_Z; j++){
          cout << "seg: " << k << ", zetaBin: " << j <<  ", zetaBinVal: " << zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j] << endl;
        }*/
      }
    }
  }

  //now map rho based on mapped z-position (binned) and segment number
  cout << "Mapping ordering parameter rho to r spatial parameter..." << endl;
  TH1 **rhoHist = new TH1*[NPOS*NCORE*NSEG*VOXEL_BINS_Z];
  TH1 **rhoHistC = new TH1*[NPOS*NCORE*NSEG*VOXEL_BINS_Z];
  TH1 **rMap = new TH1*[NPOS*NCORE*NSEG*VOXEL_BINS_Z];
  Int_t *rhoBinVal = (Int_t*)malloc((NPOS*NCORE*NSEG*(VOXEL_BINS_R + 1)*(VOXEL_BINS_Z))*sizeof(Int_t));
  memset(rhoBinVal,0,sizeof(&rhoBinVal));
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      for(int j = 0; j < VOXEL_BINS_Z; j++){
        //get the projection in rho (x) gated on a range in zeta (z)
        rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = new TH1D();
        rhoHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = new TH1D();
        Int_t lowerZetaBound = zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j];
        Int_t upperZetaBound = zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j + 1];
        Int_t zValMin = (Int_t)(j*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        Int_t zValMax = (Int_t)((j+1)*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        //cout << "seg: " << k << ", bin: " << j << ", proj bounds: " << lowerZetaBound << " " << upperZetaBound << endl;
        if((lowerZetaBound > 0)&&(lowerZetaBound < upperZetaBound)){
          rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = rhophizetaHist[l*NSEG + k]->ProjectionX("",1,N_BINS_ORDERING,lowerZetaBound,upperZetaBound);
          //cout << "proj entries: " << rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetEntries() << "orig entries: " << rhophizetaHist[l*NSEG + k]->GetEntries() << endl;
          if(rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetEntries()>0){
            //generate cumulative distributions
            rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Scale(1.0/rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Integral());
            rhoHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = rhoHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetCumulative();
            //rhoHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->SetNameTitle(Form("rhoHistCumulativePos%iSeg%ir%ito%i",l,k,rValMin,rValMax),Form("rhoHistCumulativePos%iSeg%ir%ito%i",l,k,rValMin,rValMax));
            //writeList->Add(rhoHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]);
            if(rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->Integral()>0.){
              //generate maps
              sprintf(hname,"rMapPos%iSeg%iz%ito%i",l,k,zValMin,zValMax);
              rMap[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j] = new TH1D(hname,Form("rMapPos%iSeg%iz%ito%i",l,k,zValMin,zValMax),N_BINS_ORDERING,-1.0*RHO_MAX,RHO_MAX);
              for(int m=0;m<N_BINS_ORDERING;m++){
                xVal[0] = rhoHistC[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetBinContent(m+1);
                nQuantiles = rDistHist[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->GetQuantiles(1,qVal,xVal);
                if(nQuantiles>=1){
                  //cout << "m = " << m << ", nQuantiles = " << nQuantiles << ", val = " << qVal[0] << endl; //qVal[0]=0 when there is no GEANT4 data in this bin
                  rMap[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->SetBinContent(m+1,qVal[0]);
                }
              }
              writeList->Add(rMap[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]);
              //cout << "here write rMap" << l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j << endl;
              //store bin numbers which correspond to discrete values of the r
              for(int i = 0; i <= VOXEL_BINS_R; i++){
                //printf("k j i %i %i %i\n",k,j,i);
                //printf("bins %i %i\n",k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i, k*VOXEL_BINS_Z + j);
                Double_t rBinVal = i*(MAX_VAL_R/(VOXEL_BINS_R*1.0));
                rhoBinVal[l*NSEG*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + j*(VOXEL_BINS_R+1) + i] = rMap[l*NSEG*VOXEL_BINS_Z + k*VOXEL_BINS_Z + j]->FindFirstBinAbove(rBinVal-0.01); //-0.01 to handle bin values at the end of the mapping range
                //printf("k j i %i %i %i, binVal: %i\n",k,j,i,rhoBinVal[l*NSEG*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + j*(VOXEL_BINS_R+1) + i]);
              }
              /*for(int i = 0; i <= VOXEL_BINS_R; i++){
                cout << "seg: " << k << ", rhoBin: " << i <<  ", rhoBinVal: " << rhoBinVal[l*NSEG*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + j*(VOXEL_BINS_R+1) + i] << endl;
              }*/
            }
          }
        }
      }
    }
  }

  //now map phi based on mapped z-value and radius (both binned) and segment number
  cout << "Mapping ordering parameter phi to angle spatial parameter..." << endl;
  //generate cumulative distributions
  TH1 **phiHist = new TH1*[NPOS*NCORE*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z)];
  TH1 **phiHistC = new TH1*[NPOS*NCORE*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z)];
  TH1 **angleMap = new TH1*[NPOS*NCORE*NSEG*(VOXEL_BINS_Z)*(VOXEL_BINS_R)];
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int k = 0; k < NSEG; k++){
      for(int j = 0; j < VOXEL_BINS_Z; j++){
        Int_t lowerZetaBound = zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j];
        Int_t upperZetaBound = zetaBinVal[l*NSEG*(VOXEL_BINS_Z+1) + k*(VOXEL_BINS_Z+1) + j + 1];
        Int_t zValMin = (Int_t)(j*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        Int_t zValMax = (Int_t)((j+1)*(MAX_VAL_Z/(VOXEL_BINS_Z*1.0)));
        for(int i = 0; i < VOXEL_BINS_R; i++){
          Int_t lowerRhoBound = rhoBinVal[l*NSEG*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + j*(VOXEL_BINS_R+1) + i];
          Int_t upperRhoBound = rhoBinVal[l*NSEG*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R+1)*(VOXEL_BINS_Z) + j*(VOXEL_BINS_R+1) + i + 1];
          Int_t rValMin = (Int_t)(i*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
          Int_t rValMax = (Int_t)((i+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
          //cout << "bin: " << k << " " << j << " " << i << ", proj bounds x: " << lowerBoundx << " " << upperBoundx << endl;
          //cout << "proj bounds y: " << lowerBoundy << " " << upperBoundy << endl;
          if((lowerZetaBound > 0)&&(lowerZetaBound < upperZetaBound)){
            if((lowerRhoBound > 0)&&(lowerRhoBound < upperRhoBound)){
              //get the projection in phi (y) gated on a range in rho (x) and zeta (z)
              phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = new TH1D();
              phiHistC[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = new TH1D();
              phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = rhophizetaHist[l*NSEG + k]->ProjectionY("",lowerRhoBound,upperRhoBound,lowerZetaBound,upperZetaBound);
              if(phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->GetEntries()>0){
                phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Scale(1.0/phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Integral());
                phiHistC[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->GetCumulative();
                phiHistC[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->SetNameTitle(Form("phiHistCumulativePos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax),Form("phiHistCumulativePos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax));
                //writeList->Add(phiHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]);
                //generate maps
                if(angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->Integral()>0.){
                  sprintf(hname,"angleMapPos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax);
                  angleMap[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i] = new TH1D(hname,Form("angleMapPos%iSeg%iz%ito%ir%ito%i",l,k,zValMin,zValMax,rValMin,rValMax),N_BINS_ORDERING,-1.*PHI_MAX,1.*PHI_MAX);
                  for(int m=0;m<N_BINS_ORDERING;m++){
                    xVal[0] = phiHistC[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->GetBinContent(m+1);
                    nQuantiles = angleDistHist[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->GetQuantiles(1,qVal,xVal);
                    if(nQuantiles>=1){
                      //cout << "m = " << m << ", nQuantiles = " << nQuantiles << ", xval = " << xVal[0] <<  ", qval = " << qVal[0] << endl; //qVal[0]=0 when there is no GEANT4 data in this bin
                      angleMap[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]->SetBinContent(m+1,qVal[0]);
                    }
                  }
                  //cout << "here write angleMap " << l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i << endl;
                  writeList->Add(angleMap[l*NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + k*(VOXEL_BINS_R)*(VOXEL_BINS_Z) + j*VOXEL_BINS_R + i]);
                }
              }
            }
          }
        }
      }
    }
  }

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  writeList->Write();
  myfile->Close();
  free(zetaBinVal);
  free(rhoBinVal);
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

  cout << "Input file: " << afile << endl << "Simulation data: " << simfile << endl << "Calibration file: " << calfile << endl << "Output file: " << outfile << endl;

  TParserLibrary::Get()->Load();

  bool inpList = false;
  bool simList = false;
  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/
  if(strcmp(ext,".list")==0){
    inpList = true;
    cout << "Sorting from a list of analysis trees..." << endl;
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
  }
  ext=strrchr(argv[2],'.'); /* returns a pointer to the last . to grab extention*/
  if(strcmp(ext,".list")==0){
    simList = true;
    cout << "Sorting from a list of simulation trees..." << endl;
  }else{
    cout << "Sorting from a single simulation tree..." << endl;
  }

  generate_mapping(afile, simfile, calfile, outfile, inpList, simList);
  cout << "Mapping complete!" << endl;

  return 0;
}