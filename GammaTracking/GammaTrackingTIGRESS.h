#ifndef GAMMATRACKINGTIGRESS_H
#define GAMMATRACKINGTIGRESS_H

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TTigress.h"
#include "TGriffin.h"
#include "TGRSIDetectorHit.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TPulseAnalyzer.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

#define     NPOS   16 //number of positions in the array
#define     NCORE  4  //number of cores per position
#define     NSEG   8  //number of segments per core

#define     SAMPLES          100 //number of samples in each waveform
#define     BASELINE_SAMPLES 20  //number of samples at the start of the waveform used to calculate the baseline

#define     N_BINS_ORDERING 512 //number of bins to use when discretizing ordering parameter (WARNING: memory usage scales as ^3 with this, can also overflow TH3 integer bin index with values > 1024!)
#define     RHO_MAX         5.0E4
#define     PHI_MAX         1.0
#define     ZETA_MAX        1.0

#define     MAX_VAL_R             36 //maximum r (in mm)
#define     VOXEL_BINS_R          18 //number of r bins in the map (large values can cause memory pressure when sorting the basis)
#define     MAX_VAL_ANGLE         90 //maximum angle (in deg)
#define     VOXEL_BINS_ANGLE_MAX  9 //number of angle bins in the map, at the largest r value (scales with r) (large values can cause memory pressure when sorting the basis)
#define     MAX_VAL_Z             90 //maximum z (in mm)
#define     VOXEL_BINS_Z          9 //used for binning in the basis (not needed for the map) (large values can cause memory pressure when sorting the basis)

#define     MAX_ENERGY_SINGLE_INTERACTION  1500 //maximum energy allowed for a hit to be put into the map (higher energy events are more likely to result from multiple interactions)
#define     BASIS_MAX_ENERGY        2000 //maximum energy allowed for a hit to be put into the basis (used to suppress high energy events)
#define     BASIS_MIN_ENERGY        0

//parameters defining how fine the basis grid is (number of bins in each dimension, per segment)
//bounds of the grid determined by MAX_VAL_R, MAX_VAL_ANGLE, MAX_VAL_Z
#define     COARSE_BASIS_BINFACTOR  1.0 //for the coarse basis, the number of bins in each dimension is multiplied by this factor with respect to the map
#define     FINE_BASIS_BINFACTOR    4.0 //for the fine basis, the number of bins in each dimension is multiplied by this factor with respect to the map

#define     GRID_HIT_SEG_WEIGHT     1.0
#define     GRID_NONHIT_SEG_WEIGHT  20.0

#define     SEGMENT_ENERGY_THRESHOLD       50 //threshold for a segment to be considered 'hit' in GetEnergy() units
#define     SEGMENT_ENERGY_NOHIT_THRESHOLD 20  //threshold in GetCharge() units below which a segment is considered not to be hit

#define     MAX_BASIS_SORT_CHISQ  2.0 //threshold chisq for basis sort

#define     BIG_NUMBER 1E30  //a big number
#define     BAD_RETURN -1E10 //value to be returned if ordering parameter calculation fails

typedef struct
{
  TH1 *rMap[NSEG];
  TH1 *angleMap[NSEG*VOXEL_BINS_R];
  TH1 *zMap[NSEG*(VOXEL_BINS_R)*(VOXEL_BINS_ANGLE_MAX)];
}GT_map;

typedef struct
{
  TH1D *coarseBasis[(Int_t)(VOXEL_BINS_R*COARSE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*COARSE_BASIS_BINFACTOR*VOXEL_BINS_Z*COARSE_BASIS_BINFACTOR)];
  TH1D *fineBasis[(Int_t)(VOXEL_BINS_R*FINE_BASIS_BINFACTOR*4*VOXEL_BINS_ANGLE_MAX*FINE_BASIS_BINFACTOR*VOXEL_BINS_Z*FINE_BASIS_BINFACTOR)];
  TH1I *basisHPCoarse, *basisHPFine;
}GT_basis;



//gamma tracking functions
void GT_import_map(TFile*,GT_map*);
void GT_import_basis(TFile*,TFile*,GT_basis*);
TVector3 GT_get_pos_direct(TTigressHit*,GT_map*);
TVector3 GT_get_pos_gridsearch(TTigressHit*,GT_basis*);
TVector3 GT_transform_position_to_clover(TTigressHit*,TVector3*);
TVector3 GT_transform_position_to_absolute(TTigressHit*,TVector3*);
double GT_get_doppler(double,TVector3*,TTigressHit*,TVector3*);
Int_t getMaxChargeSegHit(TTigressHit*);

#endif