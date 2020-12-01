#include "common.h" //define all global variables here!

TH1 *rMap[NSEG], *angleMap[NSEG*VOXEL_BINS_R], *zMap[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];
TH1D *basis[(Int_t)(VOXEL_BINS_R*4*VOXEL_BINS_ANGLE_MAX*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR)];

void sortData(TFile *inputfile, const char *calfile, const Double_t basisScaleFac, const bool makeFineBasis, unsigned long int *numEvtsBasis){
  
  const Int_t basisBinsR = VOXEL_BINS_R*basisScaleFac; 
  const Int_t basisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*basisScaleFac; 
  const Int_t basisBinsZ = VOXEL_BINS_Z*basisScaleFac;

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

  Int_t hit_counter = 0;
  Int_t map_hit_counter = 0;

  const Int_t one = 1;
  for (int jentry = 0; jentry < tree->GetEntries(); jentry++) {
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue; //exclude pileup
      Double_t coreCharge = tigress_hit->GetCharge();
      if((coreCharge <= 0)||(coreCharge > BASIS_MAX_ENERGY)) continue; //bad energy
      tigress_hit->SetWavefit();
      const std::vector<Short_t> *wf = tigress_hit->GetWaveform();
      if(wf->size()!=SAMPLES){
        cout << "Entry " << jentry << ", improper core waveform size (" << wf->size() << ")." << endl;
        continue;
      }
      hit_counter++;
      //cout << "Number of segments: " << tigress_hit->GetSegmentMultiplicity() << endl;
      if(tigress_hit->GetSegmentMultiplicity() == NSEG){
        //all segments have waveforms
        Int_t numSegHits = 0; //counter for the number of segments with a hit (ie. over the threshold energy)
        //check that the waveforms are the same size, and that there are no duplicate segments
        bool goodWaveforms = true;
        Int_t segsInData = 0;
        for(int i = 0; i < tigress_hit->GetSegmentMultiplicity(); i++){
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
          if((tigress_hit->GetSegmentHit(i).GetCharge() > BASIS_MAX_ENERGY)||(fabs(tigress_hit->GetSegmentHit(i).GetCharge()) > MAX_ENERGY_SINGLE_INTERACTION)){
            goodWaveforms = false;
            break;
          }
          if(tigress_hit->GetSegmentHit(i).GetCharge() > 0.2*coreCharge){
            numSegHits++;
          }
        }
        if(numSegHits != 1){
          goodWaveforms = false;
        }
        if(goodWaveforms){
          //cout << "Entry " << jentry << endl;
          bool isHit = false;

          for(int i = 0; i < tigress_hit->GetSegmentMultiplicity(); i++){

            //calculate all ordering parameters (see ordering_parameter_calc.cxx)
            double rho = calc_ordering(tigress_hit,i,jentry,0);
            if(rho == BAD_RETURN){
              continue;
            }
            double phi = calc_ordering(tigress_hit,i,jentry,1);
            if(phi == BAD_RETURN){
              continue;
            }
            double zeta = calc_ordering(tigress_hit,i,jentry,2);
            if(zeta == BAD_RETURN){
              continue;
            }
            if((fabs(rho) > RHO_MAX)||(fabs(phi) > PHI_MAX)||(fabs(zeta) > ZETA_MAX)){
              continue;
            }
   
            Int_t segNum = tigress_hit->GetSegmentHit(i).GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

            isHit = true;

            //map to spatial parameters
            if(rMap[segNum]!=NULL){
              double r = rMap[segNum]->GetBinContent(rMap[segNum]->FindBin(rho));
              Int_t rInd = (Int_t)(r*VOXEL_BINS_R/(1.0*MAX_VAL_R));
              //cout << "seg: " << segNum << ", r: " << r << ", rho: " << rho << ", ind: " << rInd << endl;
              if(rInd < VOXEL_BINS_R){
                if(angleMap[segNum*VOXEL_BINS_R + rInd]!=NULL){
                  double angle = angleMap[segNum*VOXEL_BINS_R + rInd]->GetBinContent(angleMap[segNum*VOXEL_BINS_R + rInd]->FindBin(phi));
                  if(angle>=MAX_VAL_ANGLE){
                    angle=MAX_VAL_ANGLE-0.001; //have seen rare events where angle is exactly 90 degrees
                  }
                  Int_t angleInd = (Int_t)(angle*getNumAngleBins(rInd,1.0)/MAX_VAL_ANGLE);
                  //cout << "angle: " << angle << ", phi: " << phi << ", ind: " << angleInd << endl;
                  if(zMap[segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]!=NULL){
                    double z = zMap[segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]->GetBinContent(zMap[segNum*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + rInd*VOXEL_BINS_ANGLE_MAX + angleInd]->FindBin(zeta));
                    if(z>=MAX_VAL_Z){
                      z=MAX_VAL_Z-0.001; //have seen rare events where z is exactly 90 mm
                    }
                    //cout << "segnum: " << segNum << ", r: " << r << ", angle: " << angle << ", z: " << z << endl;

                    if((30.-z) < r){
                      if((r>0.)&&(z>0.)&&(angle>0.)){
                        if(segNum<=3){
                          //r corresponds to the distance from the central contact at z=30
                          r = sqrt(r*r - (30.-z)*(30.-z));
                        }
                        if((r==r)&&(angle==angle)&&(z==z)){
                          angle += 90.*(segNum%4); //transform angle to 2pi spanned val
                          if(angle>=360.){
                            angle=360.-0.001; //have seen rare events where angle is exactly 360 degrees
                          }

                          //handle erronous mapping outside of real TIGRESS geometry due to map bin size
                          if(r<5. && z>=20.){
                            /*cout << "r: " << r << ", angle: " << angle << ", z: " << z << endl;
                            cout << "rInd: " << rInd << ", angleInd: " << angleInd << endl;
                            cout << "segment: " << segNum << ", rho: " << rho << ", phi: " << phi << ", zeta: " << zeta << endl;
                            getc(stdin);*/
                            r = 5.0;
                          }

                          //get indices for r, angle, z
                          const Int_t rBasisInd = (Int_t)(r*basisBinsR/MAX_VAL_R);
                          Int_t numAngleBinsAtR = 4*getNumAngleBins(rBasisInd,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
                          if(makeFineBasis){
                            numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
                          }
                          const Int_t angleBasisInd = (Int_t)(angle*numAngleBinsAtR/360.);
                          const Int_t zBasisInd = (Int_t)(z*basisBinsZ/MAX_VAL_Z);
                          const Int_t basisInd = rBasisInd*basisBinsAngle*basisBinsZ + angleBasisInd*basisBinsZ + zBasisInd;
                          
                          //cout << "seg: " << segNum << ", rBasisInd: " << rBasisInd << ", angleBasisInd: " << angleBasisInd << ", numAngleBinsAtR: " << numAngleBinsAtR << ", zBasisInd: " << zBasisInd << endl;
                          //cout << "num basis bins: [ " << basisBinsR << " " << basisBinsAngle << " " << basisBinsZ << " ]" << endl;

                          //save waveforms (in 'superpulse' format, core waveform followed by segment waveforms on the same histogram)

                          //first core waveform
                          Double_t core_waveform_baseline = 0.;
                          for(int k = 0; k < BASELINE_SAMPLES; k++){
                            core_waveform_baseline += wf->at(k);
                          }
                          core_waveform_baseline /= 1.0*BASELINE_SAMPLES;
                          for(int k = 0; k < SAMPLES; k++){
                            basis[basisInd]->SetBinContent(k+1,basis[basisInd]->GetBinContent(k+1) + ((wf->at(k) - core_waveform_baseline)/coreCharge));
                          }
                          //then segment waveforms
                          for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
                            //cout << "basisInd: " << basisInd << endl;

                            const Int_t segBasisInd = tigress_hit->GetSegmentHit(j).GetSegment();
                            const std::vector<Short_t> *segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                            Double_t seg_waveform_baseline = 0.;
                            for(int k = 0; k < BASELINE_SAMPLES; k++){
                              seg_waveform_baseline += segwf->at(k);
                            }
                            seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;

                            //cout << "core t0: " << waveform_t0 << ", core energy: " << coreCharge << ", segment " << segBasisInd << " baseline: " << seg_waveform_baseline << endl;
                            
                            for(int k = 0; k < SAMPLES; k++){
                              //cout << "incrementing bin " << k << "by " << ((segwf->at(k) - seg_waveform_baseline)/coreCharge) << endl;
                              basis[basisInd]->SetBinContent(k+1+(SAMPLES*segBasisInd),basis[basisInd]->GetBinContent(k+1+(SAMPLES*segBasisInd)) + ((segwf->at(k) - seg_waveform_baseline)/coreCharge) );
                            }
                          }
                          numEvtsBasis[basisInd]++;
                        }
                        //cout << "val: " << basis[basisInd]->GetBinContent(0) << ", num evts: " << numEvtsBasis[basisInd] << endl;
                      }
                    }else{
                      cout << "WARNING: lower segment invalid z value (" << z << ", r: " << r << ")" << endl;
                    }
                  }else{
                    cout << "WARNING: NULL z map bin (seg: " << segNum << ", r: " << r << ", r ind: " << rInd << ", angle: " << angle << ", angleInd: " << angleInd << ")" << endl;
                  }
                }else{
                  cout << "WARNING: NULL angle map bin (seg: " << segNum << ", r: " << r << ", r ind: " << rInd << ")" << endl;
                }
              }
            }else{
              cout << "WARNING: NULL r map bin (seg: " << segNum << ")" << endl;
            }
            
          }
          if(isHit){
            map_hit_counter++;
          }
        }
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  cout << map_hit_counter << " of " << hit_counter << " hits retained (" << 100*map_hit_counter/hit_counter << " %)." << endl;

}

//Function which generates basis waveforms from calibration data, using a map file generated by the GammaTrackingMakeMap code.
//Waveform basis data is stored in a ROOT tree.
//This assumes that the calibration data predominantly contains single hit events.
void make_waveform_basis(const char *infile, const char *mapfile, const char *calfile, const char *outfile, const bool makeFineBasis, bool inpList) {

  if(FINE_BASIS_BINFACTOR<COARSE_BASIS_BINFACTOR){
    cout << "ERROR: fine basis binning must be equal to or finer than the coarse basis binning." << endl;
    exit(-1);
  }
  if(fmod((FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR),1) > 0.){
    cout << "ERROR: the ratio FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR must be an integer." << endl;
    exit(-1);
  }

  TList * list = new TList;

  //read in histograms from map file
  TFile *inp = new TFile(mapfile,"read");
  if (!inp->IsOpen()) {
    cout << "ERROR: Could not open map file!" << endl;
    exit(-1);
  }else{
    cout << "Opened map file: " << mapfile << endl;
  }
  
  for(int k = 0; k < NSEG; k++){
    sprintf(hname,"rMapSeg%i",k);
    if((rMap[k] = (TH1*)inp->Get(hname))==NULL){
      cout << "No r coordinate map for segment " << k << endl;
    }
    for(int j = 0; j < VOXEL_BINS_R; j++){
      Int_t rValMin = (Int_t)(j*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      Int_t rValMax = (Int_t)((j+1)*(MAX_VAL_R/(VOXEL_BINS_R*1.0)));
      sprintf(hname,"angleMapSeg%ir%ito%i",k,rValMin,rValMax);
      if((angleMap[k*VOXEL_BINS_R + j] = (TH1*)inp->Get(hname))==NULL){
        cout << "No angle coordinate map for segment " << k << ", radial range " << rValMin << " to " << rValMax << endl;
      }
      Int_t numAngleBins = getNumAngleBins(j,1.0);
      for(int i = 0; i < numAngleBins; i++){
        Int_t angleValMin = (Int_t)(i*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        Int_t angleValMax = (Int_t)((i+1)*(MAX_VAL_ANGLE/(numAngleBins*1.0))); //size of phi bins depends on r
        sprintf(hname,"zMapSeg%ir%ito%iangle%ito%i",k,rValMin,rValMax,angleValMin,angleValMax);
        if((zMap[k*(VOXEL_BINS_ANGLE_MAX)*(VOXEL_BINS_R) + j*VOXEL_BINS_ANGLE_MAX + i] = (TH1*)inp->Get(hname))==NULL){
          cout << "No z coordinate map for segment " << k << ", radial range " << rValMin << " to " << rValMax << ", angle range " << angleValMin << " to " << angleValMax << endl;
        }
      }
    } 
  }

  cout << "Map file data read in." << endl;

  //setup histograms for the basis
  
  Double_t basisScaleFac;
  if(makeFineBasis){
    basisScaleFac = FINE_BASIS_BINFACTOR;
  }else{
    basisScaleFac = COARSE_BASIS_BINFACTOR;
  }
  const Int_t basisBinsR = VOXEL_BINS_R*basisScaleFac; 
  const Int_t basisBinsAngle = 4*VOXEL_BINS_ANGLE_MAX*basisScaleFac; 
  const Int_t basisBinsZ = VOXEL_BINS_Z*basisScaleFac;
  TH1I *basisHP = new TH1I("basis_hitpattern","basis_hitpattern",basisBinsR*basisBinsAngle*basisBinsZ,0,basisBinsR*basisBinsAngle*basisBinsZ);
  list->Add(basisHP);
  unsigned long int numEvtsBasis[basisBinsR*basisBinsAngle*basisBinsZ];
  memset(numEvtsBasis,0,sizeof(numEvtsBasis));
  for(int k = 0; k < basisBinsR; k++){
    Int_t numAngleBinsAtR = 4*getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
    if(makeFineBasis){
      numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
    }
    //cout << "numAngleBinsAtR: " << numAngleBinsAtR << endl;
    //printf("ind: %i, angle bins %i, actual %i, getnum: %i\n",k,basisBinsAngle,numAngleBinsAtR,getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR));
    for(int j = 0; j < numAngleBinsAtR; j++){
      for(int i = 0; i < basisBinsZ; i++){
        Int_t basisInd = k*basisBinsAngle*basisBinsZ + j*basisBinsZ + i;
        sprintf(hname,"basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/basisBinsR,(k+1)*MAX_VAL_R/basisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/basisBinsZ,(i+1)*MAX_VAL_Z/basisBinsZ);
        basis[basisInd] = new TH1D(hname,Form("basis_r%ito%i_angle%ito%i_z%ito%i",k*MAX_VAL_R/basisBinsR,(k+1)*MAX_VAL_R/basisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/basisBinsZ,(i+1)*MAX_VAL_Z/basisBinsZ),SAMPLES*(NSEG+1),0,SAMPLES*(NSEG+1));
        //list->Add(basis[basisInd]);
      }
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
      TFile * inputfile = new TFile(name, "READ");
      if (!inputfile->IsOpen()) {
        cout << "ERROR: Could not open analysis tree file (" << name << ")" << endl;
        exit(-1);
      }
      sortData(inputfile, calfile, basisScaleFac, makeFineBasis, numEvtsBasis);
      inputfile->Close();
    }
    fclose(listFile);
  }else{
    TFile * inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      cout << "ERROR: Could not open analysis tree file (" << infile << ")" << endl;
      exit(-1);
    }
    sortData(inputfile, calfile, basisScaleFac, makeFineBasis, numEvtsBasis);
    inputfile->Close();
  }
  

  for(int k = 0; k < basisBinsR; k++){
    Int_t numAngleBinsAtR = 4*getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
    if(makeFineBasis){
      numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
    } 
    for(int j = 0; j < numAngleBinsAtR; j++){
      for(int i = 0; i < basisBinsZ; i++){
        const Int_t basisInd = k*basisBinsAngle*basisBinsZ + j*basisBinsZ + i;
        uint32_t basisHPVal = 0;
        const uint32_t one = 1;
        if(numEvtsBasis[basisInd] > 0){
          for(int m = 0; m < SAMPLES*(NSEG+1); m++){
            basis[basisInd]->SetBinContent(m+1,basis[basisInd]->GetBinContent(m+1)/(1.0*numEvtsBasis[basisInd]));
          }
          for(int m = 0; m < NSEG; m++){
            //try and find the sample value near the expected maximum of the pulse (estimate it is at 0.85*SAMPLES)
            if(basis[basisInd]->GetBinContent((Int_t)((m+1.85)*(SAMPLES))) > 0.2){
              //this segment is 'hit' in this basis bin
              basisHPVal|=(one<<m);
              //cout << "hit seg " << m << " with bin content: " << basis[basisInd]->GetBinContent((Int_t)((m+1.85)*(SAMPLES))) << endl;
            }
          }
          //cout << "index: " << basisInd << ", HP val: " << basisHPVal << endl;
          //getc(stdin);
          basisHP->SetBinContent(basisInd+1,basisHPVal);
          list->Add(basis[basisInd]);
        }
      }
    }
  }

  

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  inp->Close();

}

int main(int argc, char ** argv) {

  const char *afile, *mapfile, *outfileCoarse, *outfileFine, *calfile;
  char *ext;

  if (argc < 2) {
    cout << endl << "This sortcode generates a waveform basis for gamma tracking, using the map file generated using the GammaTrackingMakeMap code." << endl << endl;
    cout << "Arguments: ./GammaTrackingMakeBasis analysis_tree_file map_file cal_file coarse_basis_output_file fine_basis_output_file" << endl << endl;
    cout << "The analysis tree (containing the calibration data used to make the basis) is a required argument.  Omitting other arguments will cause the sortcode to fall back to default values." << endl << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
    mapfile = "trackingMap.root";
	  calfile = "CalibrationFile.cal";
	  outfileCoarse = "trackingWaveformBasisCoarse.root";
    outfileFine = "trackingWaveformBasisFine.root";
  } else if (argc == 3) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = "CalibrationFile.cal";
	  outfileCoarse = "trackingWaveformBasisCoarse.root";
    outfileFine = "trackingWaveformBasisFine.root";
  } else if (argc == 4) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = argv[3];
	  outfileCoarse = "trackingWaveformBasisCoarse.root";
    outfileFine = "trackingWaveformBasisFine.root";
  } else if (argc == 5) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = argv[3];
	  outfileCoarse = argv[4];
    outfileFine = "trackingWaveformBasisFine.root";
  } else if (argc == 6) {
	  afile = argv[1];
    mapfile = argv[2];
	  calfile = argv[3];
	  outfileCoarse = argv[4];
    outfileFine = argv[5];
  } else if (argc > 6) {
	  cout << "Too many arguments." << endl;
    cout << "Arguments: ./GammaTrackingMakeBasis analysis_tree_file map_file cal_file coarse_basis_output_file fine_basis_output_file" << endl;
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

  cout << "Input file: " << afile << endl << "Map data file: " << mapfile << endl << "Calibration file: " << calfile << endl; 
  cout << "Output file (coarse basis): " << outfileCoarse << endl << "Output file (fine basis): " << outfileFine << endl;

  TParserLibrary::Get()->Load();

  //make the coarse and fine waveform bases
  ext=strrchr(argv[1],'.'); /* returns a pointer to the last . to grab extention*/
  if(strcmp(ext,".list")==0){
    cout << "Sorting from a list of analysis trees..." << endl;
    make_waveform_basis(afile, mapfile, calfile, outfileCoarse, false, true);
    make_waveform_basis(afile, mapfile, calfile, outfileFine, true, true);
  }else{
    cout << "Sorting from a single analysis tree..." << endl;
    make_waveform_basis(afile, mapfile, calfile, outfileCoarse, false, false);
    make_waveform_basis(afile, mapfile, calfile, outfileFine, true, false);
  }
  

  return 0;
}
