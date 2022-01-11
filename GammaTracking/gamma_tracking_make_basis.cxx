#include "GammaTrackingTIGRESS.h" //define all global variables here!

//forward declarations
Int_t getNumAngleBins(Int_t rInd, Double_t rScaleFac, Double_t scaleFac);
Int_t getNumAngleBins(Int_t rInd, Double_t scaleFac);
double calc_ordering(TTigressHit *, const Int_t, const Int_t);

char hname[64];
TH1D *basis[(Int_t)(NPOS*NCORE*VOXEL_BINS_R*4*VOXEL_BINS_ANGLE_MAX*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR)];
TH1I *basisHP;
unsigned long int numEvtsBasis[(Int_t)(NPOS*NCORE*VOXEL_BINS_R*4*VOXEL_BINS_ANGLE_MAX*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR*FINE_BASIS_BINFACTOR)];
GT_map *trackingMap;

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
  Long64_t nentries = AnalysisTree->GetEntries();

  TTigress * tigress = 0;
  TTigressHit * tigress_hit;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "ERROR: no TTigress branch found!" << endl;
    exit(-1);
  }

  Long64_t hit_counter = 0;
  Long64_t map_hit_counter = 0;

  for(Long64_t jentry = 0; jentry < tree->GetEntries(); jentry++){
    tree->GetEntry(jentry);
    for (int hitInd = 0; hitInd < tigress->GetMultiplicity(); hitInd++) {
      tigress_hit = tigress->GetTigressHit(hitInd);
      if(tigress_hit->GetKValue() != 700) continue; //exclude pileup
      Double_t coreCharge = tigress_hit->GetCharge();
      if((coreCharge <= BASIS_MIN_ENERGY)||(coreCharge > BASIS_MAX_ENERGY)) continue; //bad energy
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
          //cout << "Entry " << jentry << endl;
          Int_t arrayPos = tigress_hit->GetArrayNumber();
          if((arrayPos>=0)&&(arrayPos<(NPOS*NCORE))){

            bool isHit = false;

            //calculate all ordering parameters (see ordering_parameter_calc.cxx)
            double rho = calc_ordering(tigress_hit,maxChargeSeg,0);
            if(rho == BAD_RETURN){
              continue;
            }
            double phi = calc_ordering(tigress_hit,maxChargeSeg,1);
            if(phi == BAD_RETURN){
              continue;
            }
            double zeta = calc_ordering(tigress_hit,maxChargeSeg,2);
            if(zeta == BAD_RETURN){
              continue;
            }
            if((fabs(rho) > RHO_MAX)||(fabs(phi) > PHI_MAX)||(fabs(zeta) > ZETA_MAX)){
              continue;
            }
  
            Int_t segNum = tigress_hit->GetSegmentHit(maxChargeSeg).GetSegment()-1; //1-indexed from GRSIsort, convert to 0-indexed

            isHit = true;

            //map to spatial parameters
            if(trackingMap->zMap[arrayPos*NSEG + segNum]!=NULL){
              double z = trackingMap->zMap[arrayPos*NSEG + segNum]->GetBinContent(trackingMap->zMap[arrayPos*NSEG + segNum]->FindBin(zeta));
              if(z>=MAX_VAL_Z){
                z=MAX_VAL_Z-0.001; //have seen rare events where z is exactly 90 mm
              }
              Int_t zInd = (Int_t)(z*VOXEL_BINS_Z/(1.0*MAX_VAL_Z));
              //cout << "seg: " << segNum << ", z: " << z << ", zeta: " << zeta << ", ind: " << zInd << endl;
              if(zInd < VOXEL_BINS_Z){
                if(trackingMap->rMap[arrayPos*NSEG*VOXEL_BINS_Z + segNum*VOXEL_BINS_Z + zInd]!=NULL){
                  double r = trackingMap->rMap[arrayPos*NSEG*VOXEL_BINS_Z + segNum*VOXEL_BINS_Z + zInd]->GetBinContent(trackingMap->rMap[arrayPos*NSEG*VOXEL_BINS_Z + segNum*VOXEL_BINS_Z + zInd]->FindBin(rho));
                  if(r>=MAX_VAL_R){
                    r=MAX_VAL_R-0.001;
                  }
                  Int_t rInd = (Int_t)(r*VOXEL_BINS_R/MAX_VAL_R);
                  //cout << "r: " << r << ", rho: " << rho << ", ind: " << rInd << endl;
                  if(trackingMap->angleMap[arrayPos*NSEG*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + zInd*VOXEL_BINS_R + rInd]!=NULL){
                    double angle = trackingMap->angleMap[arrayPos*NSEG*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + zInd*VOXEL_BINS_R + rInd]->GetBinContent(trackingMap->angleMap[arrayPos*NSEG*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + segNum*(VOXEL_BINS_Z)*(VOXEL_BINS_R) + zInd*VOXEL_BINS_R + rInd]->FindBin(phi));
                    if(angle>=MAX_VAL_ANGLE){
                      angle=MAX_VAL_ANGLE-0.001; //have seen rare events where angle is exactly 90 degrees
                    }
                    //cout << "segnum: " << segNum << ", r: " << r << ", angle: " << angle << ", z: " << z << endl;

                    if((30.-z) < r){
                      if((r>0.)&&(z>0.)&&(angle>0.)){
                        r = getRFromREField(r,z); //transform r into cylindrical coords
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

                          //get indices for z, r, angle
                          const Int_t zBasisInd = (Int_t)(z*basisBinsZ/MAX_VAL_Z);
                          const Int_t rBasisInd = (Int_t)(r*basisBinsR/MAX_VAL_R);
                          Int_t numAngleBinsAtR = 4*getNumAngleBins(rBasisInd,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
                          if(makeFineBasis){
                            numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
                          }
                          const Int_t angleBasisInd = (Int_t)(angle*numAngleBinsAtR/360.);
                          const Int_t basisInd = arrayPos*basisBinsR*basisBinsAngle*basisBinsZ + zBasisInd*basisBinsAngle*basisBinsR + rBasisInd*basisBinsAngle + angleBasisInd;
                          
                          //cout << "seg: " << segNum << ", rBasisInd: " << rBasisInd << ", angleBasisInd: " << angleBasisInd << ", numAngleBinsAtR: " << numAngleBinsAtR << ", zBasisInd: " << zBasisInd << endl;
                          //cout << "num basis bins: [ " << basisBinsR << " " << basisBinsAngle << " " << basisBinsZ << " ]" << endl;

                          //save waveforms (in 'superpulse' format, all segment waveforms on the same histogram)
                          //then segment waveforms
                          for(int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++){
                            //cout << "basisInd: " << basisInd << endl;

                            const Int_t segBasisInd = tigress_hit->GetSegmentHit(j).GetSegment() - 1; //0-indexed segment number
                            const std::vector<Short_t> *segwf = tigress_hit->GetSegmentHit(j).GetWaveform();
                            Double_t seg_waveform_baseline = 0.;
                            for(int k = 0; k < BASELINE_SAMPLES; k++){
                              seg_waveform_baseline += segwf->at(k);
                            }
                            seg_waveform_baseline /= 1.0*BASELINE_SAMPLES;

                            //cout << "core energy: " << coreCharge << ", segment " << segBasisInd << " baseline: " << seg_waveform_baseline << endl;
                            
                            for(int k = 0; k < SAMPLES; k++){
                              if(k<numSamples){
                                //cout << "incrementing bin " << k << "by " << ((segwf->at(k) - seg_waveform_baseline)/coreCharge) << endl;
                                basis[basisInd]->SetBinContent(k+1+(SAMPLES*segBasisInd),basis[basisInd]->GetBinContent(k+1+(SAMPLES*segBasisInd)) + fabs((segwf->at(k) - seg_waveform_baseline)/coreCharge) );
                              }
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
                    cout << "WARNING: NULL angle map bin (seg: " << segNum << ", z: " << z << ", z ind: " << zInd << ", r: " << r << ", rInd: " << rInd << ")" << endl;
                  }
                }else{
                  cout << "WARNING: NULL r map bin (seg: " << segNum << ", z: " << z << ", z ind: " << zInd << ")" << endl;
                }
              }
            }else{
              cout << "WARNING: NULL z map bin (seg: " << segNum << ")" << endl;
            }
            if(isHit){
              map_hit_counter++;
            }
          }else{
            cout << "WARNING: Bad array position (" << arrayPos << ")" << endl;
          }
          
        }
      }
      
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!" << endl;
  if(hit_counter > 0){
    cout << map_hit_counter << " of " << hit_counter << " hits retained (" << 100*map_hit_counter/hit_counter << " %)." << endl;
  }else{
    cout << "ERROR: " << map_hit_counter << " of " << hit_counter << " hits retained." << endl;
    cout << "Check that the input file has valid waveform data." << endl;
    exit(0);
  }
  

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
  trackingMap = (GT_map*)malloc(sizeof(GT_map));
  TFile *mapInp = new TFile(mapfile,"read");
  GT_import_map(mapInp,trackingMap);

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
  basisHP = new TH1I("basis_hitpattern","basis_hitpattern",NPOS*NCORE*basisBinsR*basisBinsAngle*basisBinsZ,0,NPOS*NCORE*basisBinsR*basisBinsAngle*basisBinsZ);
  list->Add(basisHP);
  memset(numEvtsBasis,0,sizeof(numEvtsBasis));
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int i = 0; i < basisBinsZ; i++){
      for(int k = 0; k < basisBinsR; k++){
        Int_t numAngleBinsAtR = 4*getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
        if(makeFineBasis){
          numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
        }
        //cout << "numAngleBinsAtR: " << numAngleBinsAtR << endl;
        //printf("ind: %i, angle bins %i, actual %i, getnum: %i\n",k,basisBinsAngle,numAngleBinsAtR,getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR));
        for(int j = 0; j < numAngleBinsAtR; j++){
          const Int_t basisInd = l*basisBinsR*basisBinsAngle*basisBinsZ + i*basisBinsAngle*basisBinsR + k*basisBinsAngle + j;
          sprintf(hname,"basisPos%i_z%ito%i_r%ito%i_angle%ito%i",l,i*MAX_VAL_Z/basisBinsZ,(i+1)*MAX_VAL_Z/basisBinsZ,k*MAX_VAL_R/basisBinsR,(k+1)*MAX_VAL_R/basisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR);
          basis[basisInd] = new TH1D(hname,Form("basisPos%i_r%ito%i_angle%ito%i_z%ito%i",l,k*MAX_VAL_R/basisBinsR,(k+1)*MAX_VAL_R/basisBinsR,j*360/numAngleBinsAtR,(j+1)*360/numAngleBinsAtR,i*MAX_VAL_Z/basisBinsZ,(i+1)*MAX_VAL_Z/basisBinsZ),SAMPLES*NSEG,0,SAMPLES*NSEG);
          //list->Add(basis[basisInd]);
        }
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
  
  for(int l = 0; l < NPOS*NCORE; l++){
    for(int i = 0; i < basisBinsZ; i++){
      for(int k = 0; k < basisBinsR; k++){
        Int_t numAngleBinsAtR = 4*getNumAngleBins(k,basisScaleFac,COARSE_BASIS_BINFACTOR); //x4 since covering 2pi rather than pi/2
        if(makeFineBasis){
          numAngleBinsAtR *= FINE_BASIS_BINFACTOR/COARSE_BASIS_BINFACTOR;
        } 
        for(int j = 0; j < numAngleBinsAtR; j++){
          const Int_t basisInd = l*basisBinsR*basisBinsAngle*basisBinsZ + i*basisBinsAngle*basisBinsR + k*basisBinsAngle + j;
          uint32_t basisHPVal = 0;
          if(numEvtsBasis[basisInd] > 0){
            for(int m = 0; m < SAMPLES*NSEG; m++){
              basis[basisInd]->SetBinContent(m+1,basis[basisInd]->GetBinContent(m+1)/(1.0*numEvtsBasis[basisInd]));
            }
            for(int m = 0; m < NSEG; m++){
              //try and find the sample value near the expected maximum of the pulse (estimate it is at 0.85*SAMPLES)
              if(basis[basisInd]->GetBinContent((Int_t)((m+0.85)*(SAMPLES))) > 0.3){
                //this segment is 'hit' in this basis bin
                basisHPVal|=(1<<m);
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
  }

  cout << "Writing histograms to: " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  cout << "Histograms written, sorting complete!" << endl;
  myfile->Close();
  mapInp->Close();

}

int main(int argc, char ** argv) {

  const char *afile, *mapfile, *outfileCoarse, *outfileFine, *calfile;
  char *ext;

  if (argc < 2) {
    cout << endl << "This sortcode generates a waveform basis for gamma tracking, using the map file generated using the GammaTrackingMakeMap code." << endl << endl;
    cout << "Arguments: ./GammaTrackingMakeBasis analysis_tree_file map_file cal_file coarse_basis_output_file fine_basis_output_file" << endl << endl;
    cout << "The analysis tree (containing the calibration data used to make the basis) is a required argument.  Omitting other arguments will cause the sortcode to fall "; 
    cout << "back to default values.  The analysis tree can be a single ROOT file, or a list of ROOT files (using file extension .list) can be specified instead." << endl << endl;
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
