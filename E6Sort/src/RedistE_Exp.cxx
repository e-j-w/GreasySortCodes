#define RedistE_Exp_cxx
#include "common.cxx"
#include "RedistE_Exp.h"

using namespace std;

//make sure srand() is called somewhere
double rand_expo(const double lambda){
  double u;
  u = rand() / (RAND_MAX + 1.0);
  return -log(1 - u) / lambda;
}

//Gaussian random number generator
//using Box-Muller transform
//make sure srand() is called somewhere
double rand_gaus(double mu, double sigma) {
  static double z1;
  static int genGaus = 0;
  genGaus = !genGaus;

  if(!genGaus){
    return z1 * sigma + mu;
  }

  double u1, u2;
  do{
    u1 = (double)rand() / RAND_MAX;
    u2 = (double)rand() / RAND_MAX;
  }while (u1 <= 1e-7); // Prevent log(0)

  double z0;
  z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
  z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

  return z0 * sigma + mu;
}

void SortData(const uint8_t sp, const double enOffset, const double sigma, const double lambda){

  for(uint16_t bin=0; bin<S32K; bin++){
    if((bin % 100) == 0){cout << setiosflags(ios::fixed) << "Processing bin " << bin << "\r" << flush;}
    double binYVal = mcaIn[sp][bin];
    for(double ct = 0; ct<binYVal; ct++){
      double origE = (double)bin + (double)((rand()) / (double)(RAND_MAX)); //get a random energy between bin start and end
      double expFac = 0.0;
      if(lambda > 0.0){
        expFac = rand_expo(lambda);
      } 
      //printf("Exp: %f\n",expFac);
      //printf("Gaus: %f\n",rand_gaus(0,sigma));
      double finalE = origE + enOffset - expFac + rand_gaus(0,sigma);
      int finalBin = (int)finalE;
      if((finalBin < 0)||(finalBin > S32K)){
        finalBin = 0;
      }
      mcaOut[sp][finalBin]++;
    }
  }
  
  return;

}

int main(int argc, char **argv){

  const char *sfile;
  const char *outfile;
  double enOffset = 0.0;
  double gausWidth = 0.0;
  double expLambda = 0.0;
  int numSp = 0, spUsed = 0;
  int spToUse[MAX_SP_REDIST];
  printf("Starting RedistE_Exp\n");

  if(argc < 8){
    cout << "Redistributes energies in mca spectra using an exponential function." << endl;
    cout << "Arguments: RedistE_Exp input_dmca_file output_dmca_file en_offset gaus_width exp_lambda num_sp sp1 sp2 ..." << endl;
    cout << "  *num_sp sp1 sp2 ...* specify which spectra in the input file to redistribute." << endl;
    return 0;
  }else{
    sfile = argv[1];
    outfile = argv[2];
    enOffset = atof(argv[3]);
    gausWidth = atof(argv[4]);
    expLambda = atof(argv[5]);
    numSp = atoi(argv[6]);
    if(numSp < 1){
      cout << "ERROR: At least one spectrum must be specified." << endl;
      return 0;
    }
    if(argc < (6+numSp)){
      cout << "ERROR: Too few arguments for the number of spectra (" << numSp << ") specified." << endl;
      return 0;
    }
    if(argc < 7+MAX_SP_REDIST){
      while(spUsed < numSp){
        spToUse[spUsed] = atoi(argv[7+spUsed]);
        spUsed++;
      }
    }else{
      cout << "ERROR: Too many spectra specified (maximum " << MAX_SP_REDIST << ")." << endl;
      return 0;
    }
  }

  memset(mcaIn,0,sizeof(mcaIn)); //zero out input spectrum
  memset(mcaOut,0,sizeof(mcaOut)); //zero out output spectrum

  FILE *inp;
  if((inp = fopen(sfile, "r")) == NULL){ //open the file
    cout << "ERROR: Cannot open the input file: " << sfile << endl;
    return 0;
  }
  int num_spect=0;
  for (int i=0; i<NSPECT; i++){
    if(fread(mcaIn[i],S32K*sizeof(double),1,inp)!=1){
      num_spect=i;
      break; // dont read in more spectra than there are in the input file
    }
  }
  fclose(inp);
  cout << num_spect << " spectra read from the input file: " << sfile << endl;
  printf("Will redistribute spectra:");
  for(int i=0; i<spUsed; i++){
    printf(" %i",spToUse[i]);
  }
  printf(" with:\n  Energy offset: %f\n  Gaussian width: %f\n  Exponential tail lambda: %f\n",enOffset,gausWidth,expLambda);

  if(num_spect > 0){
    srand(92370104);
    for(uint8_t i=0; i<num_spect; i++){
      uint8_t redisted = 0;
      for(uint8_t j=0; j<numSp; j++){
        if(i == spToUse[j]){
          printf("Redistributing values in spectrum: %i\n",i);
          SortData(i,enOffset,gausWidth,expLambda);
          redisted = 1;
        }
      }
      if(redisted == 0){
        //spectrum was not redistributed, copy original to output
        memcpy(&mcaOut[i],&mcaIn[i],sizeof(mcaIn[i]));
      }
    }
    
    cout << "Writing output data to: " << outfile << endl;
    FILE *out;
    if((out = fopen(outfile, "w")) == NULL){ //open the file
      cout << "ERROR: Cannot open the output file: " << outfile << endl;
      return 0;
    }else{
      for(uint8_t i=0; i<num_spect; i++){
        fwrite(&mcaOut[i],sizeof(mcaOut[i]),1,out);
      }
      fclose(out);
    }
  }

  return 0;
}
