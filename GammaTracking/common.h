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

#define     SAMPLES 100 //number of samples in each waveform

#define     N_BINS_ORDERING 512 //number of bins to use when discretizing ordering parameter (WARNING: memory usage scales as ^3 with this, can also overflow TH3 integer bin index with values > 1024!)
#define     RHO_MAX         5E4
#define     PHI_MAX         0.5
#define     ZETA_MAX        4

#define     MAX_VAL_R       50 //maximum r (in mm)
#define     BIN_WIDTH_R     5  //in mm
#define     MAX_VAL_ANGLE   90 //maximum angle (in deg)
#define     BIN_WIDTH_ANGLE 10 //in deg
#define     MAX_VAL_Z       90 //maximum z (in mm)

//parameters defining how fine the basis grid is (number of bins in each dimension, per segment)
//bounds of the grid determined by MAX_VAL_R, MAX_VAL_ANGLE, MAX_VAL_Z
#define     BASIS_MAX_VAL_Z_FRONT  30 //maximum z for front segments (in mm)
#define     BASIS_MIN_VAL_Z_BACK   30 //minimum z for back segments (in mm)
#define     BASIS_BINS_R           10
#define     BASIS_BINS_ANGLE       36 //in the basis, the angle covers the full 2pi range
#define     BASIS_BINS_Z           10
#define     BASIS_MAX_ENERGY       2000 //maximum energy allowed for a hit to be put into the basis

#define     BAD_RETURN -1E10 //value to be returned if ordering parameter calculation fails

#define     SEGMENT_ENERGY_THRESHOLD       200 //threshold for a segment to be considered 'hit' in GetEnergy() units
#define     SEGMENT_ENERGY_NOHIT_THRESHOLD 20  //threshold in GetEnergy() units below which a segment is considered not to be hit
#define     WAVEFORM_SAMPLING_WINDOW       10  //number of waveform samples used to construct ordering parameters

#include "ordering_parameter_calc.cxx"