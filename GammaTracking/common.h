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

#define     N_BINS_ORDERING 512 //number of bins to use when discretizing ordering parameter (WARNING: memory usage scales as ^3 with this, can also overflow TH3 integer bin index with values > 1024!)
#define     RHO_MAX         1E9
#define     ZETA_MAX        1E2

#define     MAX_VAL_R       40     //maximum r (in mm)
#define     BIN_WIDTH_R     5      //in mm
#define     MAX_VAL_ANGLE   90     //maximum angle (in deg)
#define     BIN_WIDTH_ANGLE 10     //in deg

#define     BAD_RETURN -1E10 //value to be returned if ordering parameter calculation fails

#include "ordering_parameter_calc.cxx"