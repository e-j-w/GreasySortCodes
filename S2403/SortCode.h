#include "TTigress.h"
#include "TS3.h"
#include "TSRIM.h"
#include "TParserLibrary.h"
#include "TEnv.h"
#include "TCutG.h"
#include "GH1D.h"
#include "GH2D.h"

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TRandom3.h"

using namespace std;

///Change this stuff///
const double TIG_THRESH = 40.0; // Tigress energy threshold
//const double tigtigT[2] = {-175.0,175.0}; // TIGRESS - TIGRESS timing window.
const double tigtigT[2] = {-100.0,100.0}; // TIGRESS - TIGRESS timing window.
const double tigS3T[2] = {350.,450.}; // TIGRESS - S3 timing window.

//vector<string> s3_gate_names;
//vector<string> s3_gate_names = {"68Ge","196Pt","68GeUS"};
//vector<string> s3_gate_names = {"78Kr","196Pt","78KrUS"};
//vector<string> s3_gate_names = {"84Kr","196Pt","84KrUS"};
//vector<string> s3_gate_names = {"84Rb","208Pb","84KrUS"};
vector<string> s3_gate_names = {"12C","unknown"};


//68Ge
//const double beam_Z = 32.;
//const double beam_mass = 63274.6; // MeV/c^2

//78,84Kr
//const double beam_Z = 36.;
//const double beam_mass = 72582.36; // MeV/c^2
//const double beam_mass = 78163.07;

//78,80Sr
//const double beam_Z = 38.;
//const double beam_mass = 72593.4; // MeV/c^2
//const double beam_mass = 74449.2; // MeV/c^2

// 84Rb
const double beam_Z = 37;
const double beam_mass = 78166.250;

//196Pt
//const double targ_Z = 78.;
//const double targ_mass = 182540.; // MeV/c^2

//208Pb
//const double targ_Z = 82.;
//const double targ_mass = 193688.0; // MeV/c^2

//12C
const double targ_Z = 6.;
const double targ_mass = 11177.9292; // MeV/c^2

//MeV
//const double Ep = 278.8; // 68Ge bombarding energy
//const double Ep = 266.7; // 68Ge mid-target energy (196Pt)

//const double Ep = 331.5; // 78Kr bombarding energy
//const double Ep = 317.8; // 78Kr mid-target energy (196Pt)
//const double Ep = 357.0; // 84Kr bombarding energy
//const double Ep = 343.3; // 84Kr mid-target energy (196Pt)
//const double Ep = 314.0; // 80Sr "mid-target" (208Pb)
double Ep;

//const double Ep = 331.5; // 78Sr bombarding energy
//const double Ep = 340.0; // 80Sr bombarding energy

const double Ex = 0.0;

//const double target_width = 0.738; //196Pt
const double target_width = 0.882; //208Pb

//const string beam_srim_name = "ge68_in_pt196";
//const string beam_srim_name = "ge68_in_pb208";
const string beam_srim_name = "kr78_in_pt196";
//const string beam_srim_name = "kr84_in_pt196";
const string targ_srim_name = "196Pt_on_196Pt";
//const string targ_srim_name = "pb208_in_pb208";
///////////////////////


//Shouldn't need to change anything past here
const double r2d = TMath::RadToDeg();
const double d2r = TMath::DegToRad();
const double pi = TMath::Pi();
const double twopi = TMath::TwoPi();

TSRIM srimB;
TSRIM srimT;

vector<TCutG*> s3_gates;
map<string,string> dirmap;
TList* hists;
TRandom3 rand3;

TVector3 GetS3Position(int det, int ring, int sec, bool rings_facing_target = false, bool smear = false) {
  
  if(det < 1 || det > 2 || ring < 0 || ring > 23 || sec < 0 || sec > 31)
    return TVector3(0.0,0.0,0.0);
   
  TVector3 pos(1.0,1.0,0.0); //X and Y arbitrary, Z must be zero
  double rOff = 0;
  if(smear)
 	rOff = rand3.Rndm() - 0.5;
  pos.SetPerp(ring + 11.5 + rOff);

    
  double pOff = 0; 
  if(smear)
  	pOff = rand3.Rndm() * 11.5 - 11.5/2.;
 
//printf("det: %i\n",det);

  //JW: sector angles determined based on the source rotation test
  //the 9 sector offset is (I think) related to the 90 degree 
  //difference between the EMMA and Bambino chamber mounts
  //these angles are defined for detector 1 (downstream)
  //with x axis +ve in the direction of TIGRESS position 5
  pos.SetPhi((sec+9.0)*twopi/32.0);

  //S3 detectors have opposite orientation, so flip detector 2
  if(det == 2){
    pos.RotateY(pi);
  }

  //rotate both detectors based on the Bambino chamber rotation
  pos.RotateZ(-22.5*d2r + pOff*d2r);

  if(rings_facing_target)
    pos.RotateY(pi);

  pos.SetZ(33.0);

  return pos;
  
}

bool BadSeg(TTigressHit* hit) {
  
  int num = hit->GetArrayNumber();
  if(num != 5 && num != 10)
    return false;

  int seg_mult = hit->GetNSegments();
  for(int j=0;j<seg_mult;j++) {
    
    TDetectorHit seg_hit = hit->GetSegmentHit(j);
    int seg = seg_hit.GetSegment();
   
    if((seg > 4 && num == 5) || (seg == 1 && num == 10))
      return true;
    
  }

  return false;

}

bool Gate1D(double value, double min, double max) {
  
  if(min < value && value < max)
    return true;
  
  return false;
}

TCutG* LoadCut(char const *cutfile, const char* cutname) { 

  TFile cuts(cutfile,"READ");
  return (TCutG*)cuts.Get(cutname);

}

void LoadGates(vector<string> gate_names) {
  
  cout << "Loading gates..." << endl;
  s3_gates.clear();
  /*
  if(s3_gate_names.size() > 0)
    s3_gates.push_back(LoadCut("beamCut.root",s3_gate_names.at(0).c_str()));

  if(s3_gate_names.size() > 1)
    s3_gates.push_back(LoadCut("targCut.root",s3_gate_names.at(1).c_str()));

  if(s3_gate_names.size() > 2)
    s3_gates.push_back(LoadCut("beamCutUS.root",s3_gate_names.at(2).c_str()));

  cout << "Found " << s3_gates.size() << " gates" << endl;
  */

  if(s3_gate_names.size() > 0) {
    string gn = s3_gate_names.at(0);
    string fn = gn + ".root";
    
    s3_gates.push_back(LoadCut(fn.c_str(),gn.c_str()));

  }
  if(s3_gate_names.size() > 1) {
    string gn = s3_gate_names.at(1);
    string fn = gn + ".root";
    
    s3_gates.push_back(LoadCut(fn.c_str(),gn.c_str()));
  }
  if(s3_gate_names.size() > 2) {
    string gn = s3_gate_names.at(2);
    string fn = gn + ".root";
    
    s3_gates.push_back(LoadCut(fn.c_str(),gn.c_str()));
  }

 label:
  for(int i=0;i<s3_gates.size();i++) {
    
    TCutG* cut = s3_gates.at(i);
    if(!cut) {
      s3_gates.erase(s3_gates.begin()+i);
      goto label;
    }
  }
  
  cout << " found " << s3_gates.size() << " gates" << endl;
  for(auto* g : s3_gates)
    cout << " " << g->GetName() << "\n";

  return;
}

void FillHist(const char* dirname, const char* histname, int bins, double low, double high,
	      double value, double weight=1.0) {

  GH1D* hist = (GH1D*)hists->FindObject(histname);
  if(!hist) {
    
    hist = new GH1D(histname,histname,bins,low,high);
    hists->Add(hist);
    
    dirmap[histname] = dirname;
    
  }
  
  hist->Fill(value,weight);
  
  return;
}

void FillHist(const char* dirname, const char* histname, int binsX, double lowX, double highX,
	      double valueX, int binsY, double lowY, double highY, double valueY, double weight=1.0) {

  GH2D* hist = (GH2D*)hists->FindObject(histname);
  if(!hist) {
    
    hist = new GH2D(histname,histname,binsX,lowX,highX,binsY,lowY,highY);
    hists->Add(hist);
    
    dirmap[histname] = dirname;
    
  }
  
  hist->Fill(valueX,valueY,weight);

  return;
}

void WriteHists(const char* outfile = "test.root") {

  //TFile* f = new TFile(outfile,"RECREATE");
  TFile f(outfile,"RECREATE");
  for(TObject* obj : *hists) {

    GH1D* h = (GH1D*)obj;
    string dirname = dirmap[h->GetName()];
    
    if(dirname.compare("")) { 

      /*
      if(!f.FindObject(dirname.c_str()))
	      f.mkdir(dirname.c_str());

      f.cd(dirname.c_str());
      */

      int pos = dirname.find_last_of("/");
      string name = dirname;
      if(pos != string::npos)
	      name = dirname.substr(pos+1);
      
      if(!f.FindObjectAny(name.c_str()))
	      f.mkdir(dirname.c_str());
      
      f.cd(dirname.c_str());
    }
    
    h->Write();
  }
  
  f.Close();

  return;
}

//Kinematics
double Theta_CM_FP(double ThetaLAB, bool sol2 = false) { //ThetaCM from the projectile (beam) angle
  
  double tau = (beam_mass/targ_mass)/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(sin(ThetaLAB) > 1.0/tau) {

    ThetaLAB = asin(1.0/tau);
    if(ThetaLAB < 0)
      ThetaLAB += pi;  

    return asin(tau*sin(ThetaLAB)) + ThetaLAB;
  }

  if(sol2)
    return asin(tau*sin(-ThetaLAB)) + ThetaLAB + pi;
  
  return asin(tau*sin(ThetaLAB)) + ThetaLAB;
  
}

double Theta_CM_FR(double ThetaLAB, bool sol2 = false) { //ThetaCM from the recoil (target) angle
  
  double tau = 1.0/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  
  if(sin(ThetaLAB) > 1.0/tau) {

    ThetaLAB = asin(1.0/tau);
    if(ThetaLAB < 0)
      ThetaLAB += pi;

    return asin(tau*sin(ThetaLAB)) + ThetaLAB;
  }

  if(sol2)
    return -asin(tau*sin(-ThetaLAB)) - ThetaLAB;
  
  return pi - (asin(tau*sin(ThetaLAB)) + ThetaLAB);
  
}

double Theta_LAB(double thetaCM) {

  double tau = (beam_mass/targ_mass)/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = sin(thetaCM)/(cos(thetaCM) + tau);

  if(tanTheta > 0)
    return atan(tanTheta);
  
  return atan(tanTheta) + pi;

}

double Recoil_Theta_LAB(double thetaCM) {

  double tau = 1.0/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));
  double tanTheta = sin(pi - thetaCM)/(cos(pi - thetaCM) + tau);
  
  return atan(tanTheta);
  
}

double KE_LAB(double thetaCM) {

  double tau = (beam_mass/targ_mass)/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = pow(targ_mass/(beam_mass + targ_mass),2);
  double term2 = 1 + tau*tau + 2*tau*cos(thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

double Recoil_KE_LAB(double thetaCM) {

  double tau = 1.0/sqrt(1 - (Ex/Ep)*(1 + beam_mass/targ_mass));

  double term1 = beam_mass*targ_mass/pow(beam_mass + targ_mass,2);
  double term2 = 1 + tau*tau + 2*tau*cos(pi - thetaCM);
  double term3 = Ep - Ex*(1 + beam_mass/targ_mass);
  
  return term1*term2*term3;
}

void SetupKinematics(const bool proj, const TVector3& s3_pos, double& beta,  double& gam,
		     double& recon_beta, double& recon_gam, TVector3& rPos) {

  rPos.SetXYZ(0.0,0.0,1.0);
  if(proj) {
    
    double thetaCM = Theta_CM_FP(s3_pos.Theta(),false);
  
    double energy = KE_LAB(thetaCM);
    //double distance = 0.5*target_width/abs(cos(s3_pos.Theta()));
    //energy = 0.001*srimB.GetAdjustedEnergy(1000.0*energy,distance,0.01);
  
    gam = (energy/beam_mass) + 1.0;
    beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));
  
    rPos.SetTheta(Recoil_Theta_LAB(thetaCM));
    rPos.SetPhi(s3_pos.Phi() - TMath::Pi());
  
    double recon_energy = Recoil_KE_LAB(thetaCM);
    //distance = 0.5*target_width/abs(cos(rPos.Theta()));
    //recon_energy = 0.001*srimT.GetAdjustedEnergy(1000.*recon_energy,distance);
  
    recon_gam = (recon_energy)/targ_mass + 1.0;
    recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));
    
  }
  else {

    double thetaCM = Theta_CM_FR(s3_pos.Theta(),false);
    
    double energy = Recoil_KE_LAB(thetaCM);
    //double distance = 0.5*target_width/abs(cos(s3_pos.Theta()));
    //energy = 0.001*srimT.GetAdjustedEnergy(1000.0*energy,distance,0.01);
    
    gam = (energy/targ_mass) + 1.0;
    beta = TMath::Sqrt(1.0 - 1.0/(gam*gam));
    
    rPos.SetTheta(Theta_LAB(thetaCM));
    rPos.SetPhi(s3_pos.Phi() - TMath::Pi());
    
    double recon_energy = KE_LAB(thetaCM);
    //distance = 0.5*target_width/abs(cos(rPos.Theta()));
    //recon_energy = 0.001*srimB.GetAdjustedEnergy(1000.*recon_energy,distance);
    
    recon_gam = (recon_energy)/beam_mass + 1.0;
    recon_beta = TMath::Sqrt(1.0 - 1.0/(recon_gam*recon_gam));
    
  }
  
  return;
}
