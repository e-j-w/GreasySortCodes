//g++ SortCode.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData
// S1873
// SortCode.cxx
// S. A. Gillespie
// 12/12/2019
// E. J. Williams
// 6/9/2020
// D. Yates
// 6/10/2021

#define Sortcode_cxx
#include "SortCode.h"

using namespace std;

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

TApplication *theApp;


bool bDrawBraggPlots = false;
int trifBraggCounter = 0;


bool tofGate(Double_t Time)
{
  if (Time < 750 && Time > 350)
    return true; // ToF gate. Change !!!
  else
    return false;
}

bool BGtofGate(Double_t Time)
{
  if ((Time < 0 && Time > -3000) || (Time < 4000 && Time > 1000))
    return true; // Random ToF gate. Change !!!
  else
    return false;
}

bool gate1D(Double_t value, Double_t min, Double_t max)
{
  if (min < value && value < max)
    return true;
  else
    return false;
}

bool loadCutG(char const *cutfile)
{ // 2D Gate Loader. Code only uses mass gate if cut file given
  TFile *cuts = new TFile(cutfile, "READ");
  //massGate = (TCutG * ) cuts->Get("mass");
  siicGate = (TCutG *)cuts->Get("siIC");
  return true;
}

bool goodIC(double tempIC[])
{
  //  double gatemin[4] = {2400, 2400, 2400, 200}; // IC Segment Minimum energy gates. Change !!!
  //  double gatemax[4] = {3200, 3200, 3200, 3200}; // IC Segment Maximum energy gates. Change !!!
  double gatemin[4] = {400, 440, 440, 420}; // gh
  double gatemax[4] = {550, 580, 600, 580}; // gh ... for 21Ne recoils in S1873
  bool good = true;
  for (int i = 0; i < 4; i++)
  {
    if (good == false)
      break;
    if (tempIC[i] < gatemax[i] && gatemin[i] < tempIC[i])
      good = true;
    else
      good = false;
  }
  return good;
}



//vector normal to grids for trific
TVector3 *trifNormVec = new TVector3(0,-TMath::Cos(60.*TMath::DegToRad()),TMath::Sin(60.*TMath::DegToRad()));


//creating the temporary bragg peak graph to fit
TGraph createTrificBraggGraph(TTrific* event){
  TGraph tmpbragg;

  TVector3 particle = event->GetPosition();

  double ZratioX=1./trifNormVec->Dot(particle.Unit());

  // Loop over hits
  for(int h=0;h<event->GetSingMultiplicity();h++){
    TTrificHit* hit=event->GetTrificSingHit(h);
    double tE=hit->GetEnergy();
            
    // Skip "bad" hits
    if(tE<0.1)continue;
    if(hit->GetDetector()==1)continue; //skip the first grid
            
    // Calculated de/dx and "corrected Z"
    double basicZ=TTrific::fInitialSpacingCart+TTrific::fSpacingCart*(hit->GetDetector()-1);
    basicZ*=0.1;//My function estimates are in cm not mm
    double dE=tE/ZratioX;
    double x=basicZ*ZratioX;
    tmpbragg.SetPoint(tmpbragg.GetN(),x,dE);
  }           
  return tmpbragg;

}


TVector2 fitBasicBragg(TGraph *braggGraph, TF1* fitFunction, double maxEng, int maxGrid, bool bSaveGraph =  false){
  fitFunction->SetParameters(maxGrid,maxEng,1000,10000);
  fitFunction->SetRange(maxGrid-5,40);
  braggGraph->Fit(fitFunction,"QNR");

  double dedxPeak = fitFunction->GetParameter(1);
  double range = fitFunction->GetParameter(0)+sqrt(abs(dedxPeak/fitFunction->GetParameter(3)));

	
  TCanvas C1;
  gPad->Update();
  braggGraph->Fit(fitFunction,"+RQ");
  braggGraph->SetMarkerStyle(20);
  braggGraph->Draw("ap"); 

  C1.Modified();
  C1.Update();
  C1.WaitPrimitive(); 
  if (bDrawBraggPlots){
    C1.SaveAs("test.root","recreate"); 
    sleep(3);
  }
  if (bSaveGraph){
    braggGraph->SetNameTitle(Form("simpleBraggGraph_%i",trifBraggCounter),Form("Trific simple bragg graph %i;Range;dE/dx (arb.)",trifBraggCounter));
    TGraph *tempGraph = (TGraph*)braggGraph->Clone();
    tempGraph->SetNameTitle(Form("simpleBraggGraph_%i",trifBraggCounter),Form("Trific simple bragg graph %i;Range;dE/dx (arb.)",trifBraggCounter));
    trifBraggList->Add(tempGraph);
  }

  if (4000 < fitFunction->GetChisquare()){
    TGraph *temp = (TGraph*)braggGraph->Clone();
    trifBraggMulti->Add(temp);
    return TVector2(-1,-1);
  }


  return TVector2(dedxPeak,range);
}


TVector3 fitComplicatedBragg(TGraph *braggGraph, TF1 *fitFunction, TVector2 basicInterpVals, double maxEng, double maxGrid, bool bSaveGraph =  false){

  //basicInterpVals is the array from fitBasicBragg. it is of the form (dedxPeak,range)

  double sig = basicInterpVals.Y()-maxGrid;
  double xmax = maxGrid;

  fitFunction->SetParameter(0,basicInterpVals.X()*0.05);
  //fitFunction->SetParLimits(0,basicInterpVals.X()*0.04,basicInterpVals.X()*0.06);
  fitFunction->SetParameter(1,0.02);
  fitFunction->SetParameter(2,sig);
  //fitFunction->SetParLimits(2,sig*0.1,sig*2);
  fitFunction->SetParameter(3,xmax+sig*0.5);
  //fitFunction->SetParLimits(3,xmax-sig,xmax+sig*2);
  fitFunction->SetLineColor(kGreen);
  braggGraph->Fit(fitFunction,"+QNR");

  double dE0 = fitFunction->Eval(5);
  double range = fitFunction->GetParameter(3)+fitFunction->GetParameter(2)*2; //approx
  double dedxPeak = fitFunction->Eval(fitFunction->GetParameter(3)-fitFunction->GetParameter(2)*2.3); //approx

  //cout<<endl<<"\ndE0 = "<<dE0;
  //    cout<<endl<<"Z = "<<dedxPeak;
  //   cout<<endl<<"range = "<<range<<endl;
	
  TCanvas C1;
  gPad->Update();
  braggGraph->Fit(fitFunction,"+RQ");
  braggGraph->SetMarkerStyle(20);
  braggGraph->Draw("ap");
 
  C1.Modified();
  C1.Update();
  C1.WaitPrimitive();
  if (bDrawBraggPlots){
    C1.SaveAs("test2.root","recreate"); 
    sleep(3);
  }

  if (bSaveGraph){
    braggGraph->SetNameTitle(Form("complexBraggGraph_%i",trifBraggCounter),Form("Trific complex bragg graph %i;Range;dE/dx (arb.)",trifBraggCounter++));
    TGraph *tempGraph = (TGraph*)braggGraph->Clone();
    tempGraph->SetNameTitle(Form("complexBraggGraph_%i",trifBraggCounter),Form("Trific complex bragg graph %i;Range;dE/dx (arb.)",trifBraggCounter++));
    trifBraggList->Add(tempGraph);
  }


  return TVector3(dedxPeak,range,dE0);

}




double pi = TMath::Pi();
double tigtigT[2] = {-100, 100}; // TIGRESS - TIGRESS Timing. Change !!!

double beta;
double thetalab;
double exc;
double ekin;
double recoiltheta;
double thetacm;
double rekin;

bool suppTig = false;
bool suppAdd = false;

bool bNoTigFilling = false; //debugging which skips tigress event filling

void SortCode::SortData(char const *afile, char const *calfile, char const *outfile, char const *target = "NULL")
{
  Initialise();

  TFile *analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen())
    {
      printf("Opening file %s failed, aborting\n", afile);
      return;
    }

  printf("File %s opened\n", afile);
  TChain *AnalysisTree = (TChain *)analysisfile->Get("AnalysisTree");
  long analentries = AnalysisTree->GetEntries();
  const char *testval = "NULL";

  // Checks for branches and sets pointers
  
  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress"))
    {
      AnalysisTree->SetBranchAddress("TTigress", &tigress);
    }
  else
    {
      cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
    }

  TSharc *sharc = 0;
  if (AnalysisTree->FindBranch("TSharc"))
    {
      AnalysisTree->SetBranchAddress("TSharc", &sharc);
    }
  else
    {
      cout << "Branch 'TSharc' not found! TSharc variable is NULL pointer" << endl;
    }

  TTrific *trif = 0;
  if (AnalysisTree->FindBranch("TTrific"))
    {
      AnalysisTree->SetBranchAddress("TTrific", &trif);
    }
  else
    {
      cout << "Branch 'TTrific' not found! TTrific variable is NULL pointer" << endl;
    }


  TSRIM *srim_93SrD = new TSRIM;
  srim_93SrD->ReadEnergyLossFile("93Sr_CD2.txt");
  TSRIM *srim_93SrH = new TSRIM;
  srim_93SrH->ReadEnergyLossFile("93Sr_CH2.txt");

  TSRIM *srim_86KrD = new TSRIM;
  srim_86KrD->ReadEnergyLossFile("86Kr_CD2.txt");
  TSRIM *srim_86KrH = new TSRIM;
  srim_86KrH->ReadEnergyLossFile("86Kr_CH2.txt");

  double EBeam = 8.0 * 92.914024; // 93Sr beam
  printf("Beam energy: %f MeV\n", EBeam);

  // Change below with actual target thicknesses
  if (strcmp(target, "two") == 0){
    printf("Using target %s - 449 \u03BCm CD2 \n", target); // Assume reaction happens in middle
//    EBeam = srim_86KrD->GetAdjustedEnergy(EBeam * 1000, 449./2., 0.001) / 1000.;
    EBeam = srim_93SrD->GetAdjustedEnergy(EBeam * 1000, 2.151/2., 0.001) / 1000;               // GetAdjustedEnergy(energy, target thickness um,step size)

  }else if (strcmp(target, "three") == 0){
    printf("Using target %s - 369 \u03BCm CD2  \n", target);
//    EBeam = srim_86KrD->GetAdjustedEnergy(EBeam * 1000, 369./2., 0.001) / 1000.;
    EBeam = srim_93SrD->GetAdjustedEnergy(EBeam * 1000, 2.656/2., 0.001) / 1000; // GetAdjustedEnergy(energy, target thickness um,step size)
 
  }else if (strcmp(target, "four") == 0){
    printf("Using target %s - 164 \u03BCm CD2  \n", target);
//    EBeam = srim_86KrD->GetAdjustedEnergy(EBeam * 1000, 164./2., 0.001) / 1000,;
     EBeam = srim_93SrD->GetAdjustedEnergy(EBeam * 1000, 164./2., 0.001) / 1000; // GetAdjustedEnergy(energy, target thickness um,step size)

  }else{
    printf("No target specified, not adjusting beam energy.\n");
  }
  printf("Adjusted Beam energy: %f MeV\n", EBeam);

  //this sets TTrific::fTargetToWindowCart distance for the TRIFIC-SHARC setup to the correct value. This is needed for good trific position reconstruction
  if (trif) trif->SetSharc(); //use this for the experiment!
  //if (trif) trif->SetSharc(0); //this is used for testing on runs 51994. Should be degrader (aka FC) to window distance. No idea what it is
  

  TReaction *sr93_rxn = new TReaction("sr93", "d", "p", "sr94", EBeam, 0, true);

  //Defining Pointers
  TTigressHit *tig_hit, *add_hit, *add_hit2;
  TSharcHit *barrel_hit, *quad_hit, *sharc_hit;
  TTrificHit *trif_sing_hit, *trif_x_hit, *trif_y_hit;


  TVector3 s3pos, recoil_vec;


  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  /*cout << "TIGRESS positions: " << endl;
    for(int det=1;det<17;det++){
    for(int cry=0;cry<4;cry++){
    TVector3 pos = tigress->GetPosition(det, cry, 0, 110., false);
    cout << "det: " << det << ", cry: " << cry << ", position: [ " << pos.x() << " " << pos.y() << " " << pos.z() << " ]" << endl;
    }
    }*/
 
  int	tig_emma_counter = 0; 
  int	tig_counter = 0;
  int emma_counter = 0;

  int	counter_ssb = 0;
  int	counter_na22 = 0;

  double_t betaDoppler = 0.1303; //beta for Doppler correction (93Sr rough)



  /*
    Things to add:
    TRIFIC:
    interpolated range vs interpolated bragg peak
    interpolated range vs interpolated de/dx @ window
    interpolated bragg peak vs interpolated de/dx at window
    de/dx vs Z
    SHARC:
    3d hitmap
    de/dx vs total E (correcting for dE detector thickness as a function of angle)
    excitation energy vs ??	
    projection Y for different sections of SHARC
    Figure 4 in proposal, but with energies on y (gammas gated on 2+ -> g.s.)


  */


  //defining functions to fit the trific spectrum with
  TF1 *peak = new TF1("peak","[1]-(x<[0])*abs([2])*pow(x-[0],2)-(x>=[0])*abs([3])*pow(x-[0],2)"); //fitting bragg peak
  TF1 *newbragg = new TF1("newbragg","[0]*exp(([1]*[1]*[2]*[2]*0.5)+([1]*(x-[3])))*TMath::Erfc(([1]*[2]/sqrt(2))-(([3]-x)/([2]*sqrt(2))))",0,50);//good approx fror bragg curve

  peak->SetParameters(320,20000,100,1000);//initial guess for peak parameters. We don't need to set this every time, since after the first instance of it converging, the values will be close enough for the next event
    

  //TFile* cutFile = new TFile("new_cut.root");
  TFile* cutFile = new TFile("sharc_E_theta_proton_cut.root");
  TCutG* cut = (TCutG*)cutFile->Get("_cut0");

  TFile* cutFileT = new TFile("TrifCut1.root");
  TCutG* cutT = (TCutG*)cutFileT->Get("_cut0");
  
  //TFile* cutFile_tmp = new TFile("test_cut_tmp.root");
  //TCutG* cut_tmp = (TCutG*)cutFile_tmp->Get("_cut0");

  printf("\nSorting analysis events...\n");
  //std::cout << "\bWarning: only sorting 1% of the data";
  for (int jentry = 0; jentry < analentries; jentry++)
  
    {
      //tigress = NULL;
      //trif = NULL;
      AnalysisTree->GetEntry(jentry);

      if(tigress)
	if(tigress->GetMultiplicity() > 0)
	  tig_counter++;


      if (tigress)
	{	
	  for (int t = 0; t < tigress->GetMultiplicity(); t++){
	    if (bNoTigFilling) continue;
	    tig_hit = tigress->GetTigressHit(t);
	
	    if(tig_hit != NULL){
	      suppTig = tig_hit->BGOFired();
	      if (!suppTig && tig_hit->GetEnergy() > 15){
		tigE->Fill(tig_hit->GetEnergy());
		tigE_ANum->Fill(tig_hit->GetArrayNumber(), tig_hit->GetEnergy());
		tigE_theta->Fill(tig_hit->GetEnergy(),TMath::RadToDeg()*(tig_hit->GetPosition().Theta()));
						
		tigDopp->Fill(tig_hit->GetDoppler(0.11));
		tigDopp_ANum->Fill(tig_hit->GetArrayNumber(), tig_hit->GetDoppler(0.11));
	      }
	    }
        
	  }

	  for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
	    if (bNoTigFilling) continue;
	    add_hit = tigress->GetAddbackHit(t);
	    suppAdd = add_hit->BGOFired();
	    if (!suppAdd && add_hit->GetEnergy() > 15){
	      addE->Fill(add_hit->GetEnergy());
	      addDoppRaw->Fill(add_hit->GetDoppler(betaDoppler));
	      addE_ANum->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
	      num_addr->Fill(add_hit->GetArrayNumber(),add_hit->GetAddress());

	      //TIGRESS-TIGRESS addback
	      for (int t2 = t+1; t2 < tigress->GetAddbackMultiplicity(); t2++){
		add_hit2 = tigress->GetAddbackHit(t2);
		suppAdd = add_hit2->BGOFired();
		if (!suppAdd && add_hit2->GetEnergy() > 15){
		  addT_addT->Fill(add_hit->GetTime() - add_hit2->GetTime());
		  addE_addE->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
		  addE_addE->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
		  if (gate1D((add_hit->GetTime() - add_hit2->GetTime()), tigtigT[0], tigtigT[1])){
		    addE_addE_tg->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
		    addE_addE_tg->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
		  }
		}
	      }
          
	    }
	  }

	}
      if (sharc){

	for (int i = 0; i < sharc->GetMultiplicity(); i++){
	  sharc_hit = sharc->GetSharcHit(i);
	  double eSharc = sharc_hit->GetDeltaE();
	  int detSharc = sharc_hit->GetDetector();
	  double frontSharc = sharc_hit->GetFrontStrip(); //revist front vs back strip
	  double backSharc = sharc_hit->GetBackStrip();
	  TVector3 posSharc = sharc_hit->GetPosition();


	  sharcE->Fill(eSharc);
	  sharcE_DetNum->Fill(detSharc,eSharc);
	  sharcE_FB->Fill(eSharc,sharc_hit->GetDeltaBackE());
	  if(abs(eSharc - sharc_hit->GetDeltaBackE()) < 200)
		  sharcT_FB->Fill(sharc_hit->GetDeltaT() - sharc_hit->GetDeltaBackT());

	  if(abs(eSharc - sharc_hit->GetDeltaBackE()) > 300 || abs(sharc_hit->GetDeltaT() - sharc_hit->GetDeltaBackT()) > 200)
		continue;

	  sharcE_front->Fill(100*(detSharc)+frontSharc,eSharc); //this will offset them so that detector one is channels 0-23, det 2 is 100-123, etc.  multiply by 100 starts each new detector at a hunded
	  //sharcE_back->Fill(100*(detSharc)+backSharc,eSharc); //same logic as for fronts
	  sharcE_back->Fill(100*(detSharc)+backSharc,sharc_hit->GetDeltaBackE());
	  //sharcE_theta->Fill(TMath::RadToDeg()*posSharc.Theta(),eSharc);
	  sharcE_theta->Fill(TMath::RadToDeg()*posSharc.Theta(),sharc_hit->GetEnergy());
	  sharcEB_theta->Fill(TMath::RadToDeg()*posSharc.Theta(),sharc_hit->GetDeltaBackE());
	  sharcE_phi->Fill(TMath::RadToDeg()*posSharc.Phi(),eSharc);
	  sharcE_theta_array[detSharc]->Fill(TMath::RadToDeg()*posSharc.Theta(),eSharc);
       
	  sharc_ffE->Fill(sharc_hit->GetFrontAddress(),sharc_hit->GetFrontCharge());
	  sharc_bfE->Fill(sharc_hit->GetBackAddress(),sharc_hit->GetBackCharge());	
 
	  //3d hitmap of all sharc events
	  sharc_3D_hitmap->Fill(posSharc.X(),posSharc.Y(),posSharc.Z());

	  //filling hitmap projections for all three projection dimenionsions (XY, XZ, YZ) both overall and detector by detector
	  sharc_XY_array[detSharc]->Fill(posSharc.X(),posSharc.Y());
	  sharc_XZ_array[detSharc]->Fill(posSharc.X(),posSharc.Z());
	  sharc_YZ_array[detSharc]->Fill(posSharc.Y(),posSharc.Z());
	  sharc_XY_projection->Fill(posSharc.X(),posSharc.Y());
	  sharc_XZ_projection->Fill(posSharc.X(),posSharc.Z());
	  sharc_YZ_projection->Fill(posSharc.Y(),posSharc.Z());
	   
	  //std::cout << "\nTrying to fill the dE E curve now";fflush(stdout);
	  //plotting pad number vs dssd number
	  sharc_pad_vs_dssd->Fill(detSharc, sharc_hit->GetPad().GetDetector());


	  
	   

	  //plot the de vs E curve
	  //do this only if there is a pad hit
	  if (0 != sharc_hit->GetPad().GetDetector()){
	    sharc_dE_E->Fill(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE());
	    sharc_dE_E_array[detSharc]->Fill(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE()); //
	     
	    //debugging plots for how the GetEnergy stuff work
	    sharcDebug_dE_E->Fill(sharc_hit->GetEnergy(),sharc_hit->GetDeltaE());
	    sharcDebug_dE_PadE->Fill(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE());
	    sharcDebug_E_PadE->Fill(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetEnergy());
	    sharcDebug_dEPadE_E->Fill(sharc_hit->GetEnergy(),sharc_hit->GetPad().GetEnergy()+sharc_hit->GetDeltaE());
	  } 
	   
	
	  //making the negative and positive projections, since the boxes overlap with each other in projections
	  if (posSharc.X() < 0){
	    sharc_YZ_projection_neg->Fill(posSharc.Y(),posSharc.Z());
	  }else {
	    sharc_YZ_projection_pos->Fill(posSharc.Y(),posSharc.Z());
	  }
	  if (posSharc.Y() < 0){
	    sharc_XZ_projection_neg->Fill(posSharc.X(),posSharc.Z());
	  }else {
	    sharc_XZ_projection_pos->Fill(posSharc.X(),posSharc.Z());
	  }
	  if (posSharc.Z() < 0){
	    sharc_XY_projection_neg->Fill(posSharc.X(),posSharc.Y());
	  }else {
	    sharc_XY_projection_pos->Fill(posSharc.X(),posSharc.Y());
	  }



//	  double excEng = sr93_rxn->GetExcEnergy(eSharc * 1e-3, sharc_hit->GetTheta(sharc->GetXOffset(),sharc->GetYOffset(),sharc->GetZOffset()),2); // 1e-3 converts from keV to MeV. integer at the end does nothing?
	  double excEng = sr93_rxn->GetExcEnergy(eSharc * 1e-3, sharc_hit->GetPosition().Theta(),2); // 1e-3 converts from keV to MeV. integer at the end does nothing?
          if(sharc_hit->GetPosition().Theta()*TMath::RadToDeg() > 90)
		  sharc_excE->Fill(excEng);
	  sharc_excE_theta->Fill(sharc_hit->GetPosition().Theta()*TMath::RadToDeg(),excEng);

	  //sharc-tigress coincidences
	  if (tigress){
	    for (int tigH = 0; tigH < tigress->GetAddbackMultiplicity(); tigH++){
	      add_hit = tigress->GetAddbackHit(tigH);
	      /*if(gate1D(add_hit->GetEnergy(), 40, 80)){
		tigSharcT->Fill(add_hit->GetTime() - sharc_hit->GetTime());
		}*/
	      tigSharcT->Fill(add_hit->GetTime() - sharc_hit->GetTime());
	      tigSharcTE->Fill(add_hit->GetTime() - sharc_hit->GetTime(),add_hit->GetEnergy());
	      tigSharcT_front->Fill(100*(detSharc)+frontSharc,add_hit->GetTime() - sharc_hit->GetTime());
	      tigSharcT_back->Fill(100*(detSharc)+backSharc,add_hit->GetTime() - sharc_hit->GetTime());
	      tigSharcT_phi->Fill(TMath::RadToDeg()*posSharc.Phi(),add_hit->GetTime() - sharc_hit->GetTime());
	      if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(), 0, 400)){
		//time gated SHARC-TIGRESS
		sharcEtigE->Fill(eSharc,add_hit->GetEnergy());
		tigE_sharc_gated->Fill(add_hit->GetEnergy());
	      }
	      if(gate1D(add_hit->GetEnergy(), 56, 63)){
		sharcE_tig_gated->Fill(eSharc);
	      }

	      if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),0, 400)) {
		//gates on sharc transfer channels:
		if (0 != sharc_hit->GetPad().GetDetector() && detSharc > 4 && detSharc < 9){
		  if ( pCut_array[detSharc]->IsInside(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE())){
		    tigE_sharc_p_gated->Fill(add_hit->GetDoppler(0.115));
		  }
		  else if (dCut_array[detSharc]->IsInside(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE())){
		    tigE_sharc_d_gated->Fill(add_hit->GetDoppler(0.115));
		  }
		  else if (tCut_array[detSharc]->IsInside(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE())){
		    tigE_sharc_t_gated->Fill(add_hit->GetDoppler(0.115));
		  }
		  /*if (aCut_array[detSharc]->IsInside(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE())){
		    tigE_sharc_a_gated->Fill(add_hit->GetEnergy());
		    }*/
		  if (cCut_array[detSharc]->IsInside(sharc_hit->GetPad().GetEnergy(),sharc_hit->GetDeltaE())){
		    tigE_sharc_c_gated->Fill(add_hit->GetDoppler(0.115));
		  }


		}

	      }

            
	      //if(cut_tmp->IsInside(TMath::RadToDeg()*posSharc.Theta(),eSharc)) {		    
	      //		tig_sharc_cut_dopp->Fill(add_hit->GetDoppler(0.11));
	      //	}

            
	      //if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),0, 400)) {
	      if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),100, 270)) {
		//if(eSharc > 500. || sharc_hit->GetBack().GetEnergy() > 500. || detSharc > 12) {
		if(cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),eSharc) || cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),sharc_hit->GetDeltaBackE())){		    
		  tigDop_gated->Fill(add_hit->GetDoppler(0.115));

		  for(int i=0;i<20;i++) {
		    double beta = (i+1)*0.01;
		    tigDop_gated_scan->Fill(i,add_hit->GetDoppler(beta));
		  }	
		}
	      }
	      //if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),0, 400)) {
	      if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),-170, 0)) {
		//if(eSharc > 500. || sharc_hit->GetBack().GetEnergy() > 500. || detSharc > 12) {
		if(cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),eSharc) || cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),sharc_hit->GetDeltaBackE())){		    
		  tigDop_gated_bkgn_low->Fill(add_hit->GetDoppler(0.115));

		  /*for(int i=0;i<20;i++) {
		    double beta = (i+1)*0.01;
		    tigDop_gated_scan->Fill(i,add_hit->GetDoppler(beta));
		  }*/	
		}
	      }
	      //if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),0, 400)) {
	      if(gate1D(add_hit->GetTime() - sharc_hit->GetTime(),330, 500)) {
		//if(eSharc > 500. || sharc_hit->GetBack().GetEnergy() > 500. || detSharc > 12) {
		if(cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),eSharc) || cut->IsInside(TMath::RadToDeg()*posSharc.Theta(),sharc_hit->GetDeltaBackE())){		    
		  tigDop_gated_bkgn_high->Fill(add_hit->GetDoppler(0.115));

		  /*for(int i=0;i<20;i++) {
		    double beta = (i+1)*0.01;
		    tigDop_gated_scan->Fill(i,add_hit->GetDoppler(beta));
		  }*/	
		}
	      }

	
	    }
	  }

	}

	 

      }

      if (trif){

	//bool for keeping track of if we've filled the single event hitmap or not
	bool bFilledSingle = false;

	double trifBraggE = 0; //for tracking the event with max energy deposit in a grid
	int trifBraggEDet = 0;
	//bool to keep track if we've had a good reconstruction or not
	bool bTrifGoodRecon = false;
	//Z->R ratio for better dE/dR reconstruction
	double ZratioX;

	//TGraph for temporary dE/dR data points
	TGraph tempBragg;

	double	tmpE = 0;
	double 	tmpdE = 0;
	double	vetoE = 0;
	
	for (int i = 0; i < trif->GetSingMultiplicity(); i++){
	  trif_sing_hit = trif->GetTrificSingHit(i);

	  double eng = trif_sing_hit->GetEnergy();
	  int det = trif_sing_hit->GetDetector();

	  if (99 != det) {
	    if(tigress) {
	      for (int tigH = 0; tigH < tigress->GetAddbackMultiplicity(); tigH++){
		add_hit = tigress->GetAddbackHit(tigH);
		TTTime->Fill(add_hit->GetTime() - trif_sing_hit->GetTime());	

		if(gate1D(add_hit->GetTime() - trif_sing_hit->GetTime(),380,430)) {
		  if(cutT->IsInside(det,eng)) {
		    trif_cut1_TigE->Fill(add_hit->GetDoppler(0.115));

		    for(int i=0;i<20;i++) {
		      
		      double beta = (i+1)*0.01;
		      trif_tig_gated_scan->Fill(i,add_hit->GetDoppler(beta));
		    }
		  }
		  
		}
				
	      }
	  
	    }
	  }

	  if(det > 17)
		vetoE += trif_sing_hit->GetEnergy();
	  if(det > 2 && det < 20)
		tmpE+=trif_sing_hit->GetEnergy();
	  if(det > 2 && det < 10)
		tmpdE+=trif_sing_hit->GetEnergy();

	  trifE->Fill(trif_sing_hit->GetEnergy());
	  trifChan->Fill(det);
	  trifE_chan->Fill(det,eng);

	  if(cutT->IsInside(det,eng)) {
	    trifE_chan_cut->Fill(det,eng);
	  }

	  if(sharc){
	    for (int sh = 0; sh < sharc->GetMultiplicity(); sh++){
	      sharc_hit = sharc->GetSharcHit(sh);
	      //cout << "tdiff: " << trif_sing_hit->GetTime() - sharc_hit->GetTime() << endl;
	      trificSharcT->Fill(trif_sing_hit->GetTime() - sharc_hit->GetTime());
	    }
	  }

	  if (eng > trifBraggE){//this will locate the grid with the max energy deposit in it for each event. We'll fill the histogram after the trific event loop

	    trifBraggE = eng;
	    trifBraggEDet = det; 
	  }

	  TVector3 trifPosition = trif->GetPosition(det); //get the position at the current grid
	   
	  //check if we had bad reconstruction. We will document the number of good vs bad reconstructions
	  if (0 == i){ //only do this once per event
		
	    if (TVector3(-100,-100,-100) == trifPosition){ //this is a bad reconstruction vector. It will appear as (-100,-100,-100) regardless of which detector you are calling GetPosition() at.
			
	      trif_XYStatus->Fill(3);	//bin 3 in this histogram is the bad reconstruction bin
			
	      //break;
	      continue;
	    } else {
	      bTrifGoodRecon = true;
	      trif_XYStatus->Fill(1); //bin 1 is the good recon bin
	      trif_GOOD_hitmap->Fill(trifPosition.X(),trifPosition.Y());
	    }
	  }

	  if (!bTrifGoodRecon) continue; // we don't want to do anything more if we had a bad reconstruction
	  if (!bFilledSingle){ //fill this single hitmap only once per event
	    trif_single_hitmap->Fill(trifPosition.X(),trifPosition.Y());
	    bFilledSingle = true;
	  }

	    
	    
	  //since we have a good reconstruction event, calculate the Z->R ratio
	  ZratioX = 1./trifNormVec->Dot(trifPosition.Unit());
	  //store this data point in the temp bragg TGraph
	  tempBragg.SetPoint(tempBragg.GetN(),det,eng);	
		
	  if (99 == det) continue; //this is because the single grid readouts for the multi-grid locations are set to channel 99 b/c they shouldn't have any events. But if they do, they will show up on channel 99, which will seg-fault the array.
	  //if (-60 > trifPosition.X()) std::cout << "\nTrif position is " << trifPosition.Print();

	  //if (0 > trifPosition.Z()) continue;
	  trif_hitmap_array[det]->Fill(trifPosition.X(),trifPosition.Y());

	  trif_hitmap->Fill(trifPosition.X(),trifPosition.Y()); //this is the overall hitmap of every event at every grid
	  trif_3D_hitmap->Fill(trifPosition.X(),trifPosition.Y(),trifPosition.Z());
	  trif_proj_YZ->Fill(trifPosition.Z(),trifPosition.Y()); //yes, it should be (Z,Y)
	  trif_proj_XZ->Fill(trifPosition.Z(),trifPosition.X()); //yes, it should be (Z,X)

	  if ( -40 > trifPosition.X()) trif_BAD_hitmap->Fill(trifPosition.X(),trifPosition.Y());
	   
 
	  if(sharc){
	    for (int sh = 0; sh < sharc->GetMultiplicity(); sh++){
	      sharc_hit = sharc->GetSharcHit(sh);
	      //cout << "tdiff: " << trif_sing_hit->GetTime() - sharc_hit->GetTime() << endl;
	      trificSharcT_GOOD->Fill(trif_sing_hit->GetTime() - sharc_hit->GetTime());
	    }
	  }
	  
	  
	}

	if(tmpE >0 && tmpdE>0 && tmpE != tmpdE && vetoE < 2000)
		trifEdE->Fill(tmpdE,tmpE);	


	if (tigress){
		//tig gated on trific PID plots
		if (sr93PIDCut->IsInside(tmpdE,tmpE)){
			trifEdE_93Sr_accepted->Fill(tmpdE,tmpE);
			for (int tigH = 0; tigH < tigress->GetAddbackMultiplicity(); tigH++){
				add_hit = tigress->GetTigressHit(tigH);
				trif_cut_93SrPID_tigE->Fill(add_hit->GetEnergy());
				trif_cut_93SrPID_DoppE->Fill(add_hit->GetDoppler(betaDoppler));
			}
			
		}
		if (rb93PIDCut->IsInside(tmpdE,tmpE)){
			tridEdE_93Rb_accepted->Fill(tmpdE,tmpE);
			for (int tigH = 0; tigH < tigress->GetAddbackMultiplicity(); tigH++){
				add_hit = tigress->GetTigressHit(tigH);
				trif_cut_93RbPID_tigE->Fill(add_hit->GetEnergy());
				trif_cut_93RbPID_DoppE->Fill(add_hit->GetDoppler(betaDoppler));
			}
			
		}
	}




	//we only want to do these full-event analysis if we had a good event reconstruction:
	//this takes a while, so only do the analysis once every 100 events, and only save the graphs once every 1000 events (up to 100 graphs)	
	//16 Nov. I (DY) took this out b/c we're running a sharc front-back filter and a 10x (or 100x) attenuator into trific. So this already limits the trific rates a lot

	//17 Nov. DY took out to speed up sort code
	TVector2 simpleBraggVals;
	TVector3 complexBraggVals;

	if (bTrifGoodRecon){// && 0 == jentry%100 && 100 > trifBraggCounter){

	  TGraph tempBragg = createTrificBraggGraph(trif);
	  if (0 == tempBragg.GetN()) continue; //skip the fitting if the graph is empty
	  simpleBraggVals = fitBasicBragg(&tempBragg, peak, trifBraggE, trifBraggEDet);//, true); //simpleBraggVals is the array from fitBasicBragg. it is of the form (dedxPeak,range)
	  //complexBraggVals = fitComplicatedBragg(&tempBragg, newbragg, simpleBraggVals, trifBraggE, trifBraggEDet);//, true); // complexBraggVals is TVector3(dedxPeak,range,dE0)

	}/*else if (bTrifGoodRecon){// && 0 == jentry%100){ //DY Same comment as above
		
		TGraph tempBragg = createTrificBraggGraph(trif);
		simpleBraggVals = fitBasicBragg(&tempBragg, peak, trifBraggE, trifBraggEDet); //simpleBraggVals is the array from fitBasicBragg. it is of the form (dedxPeak,range)
		complexBraggVals = fitComplicatedBragg(&tempBragg, newbragg, simpleBraggVals, trifBraggE, trifBraggEDet); // complexBraggVals is TVector3(dedxPeak,range,dE0)

		}*/



	trifdEdx_Range->Fill(trif->GetRange(),trifBraggE);

	if (-1 != simpleBraggVals.X()){	

	  trifdEdx_Range_interp_simple->Fill(simpleBraggVals.Y(),simpleBraggVals.X());
	  //trifdEdx_Range_interp_complex->Fill(complexBraggVals.Y(),complexBraggVals.X());
	}


	//gate on the cuts and make a few basic tig and sharc histos.

	/*if (trifCut1->IsInside(simpleBraggVals.Y(),simpleBraggVals.X()){
	  if (tigress){
	  for (int tigH = 0; tigH < tigress->GetAddbackMultiplicity(); tigH++){
	  add_hit = tigress->GetAddbackHit(tigH);

				
	  }
	  }
	  }*/

      }

      if (jentry % 10000 == 0)
	cout << setiosflags(ios::fixed) << "Entry " 
	     << jentry << " of " << analentries << ", " 
	     << 100 * jentry / analentries << "% complete" << "\r" 
	     << flush;

    } // analysis tree
  
  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;


  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();


  TDirectory *tgatedir = myfile->mkdir("Time_Gated_TIGRESS");
  tgatedir->cd();
  tgList->Write();
  myfile->cd();

  TDirectory *rtgatedir = myfile->mkdir("Random_Time_Gated_TIGRESS");
  rtgatedir->cd();
  bgList->Write();
  myfile->cd();

  TDirectory *pidgatedir = myfile->mkdir("PID_Gated_TIGRESS");
  pidgatedir->cd();
  PIDgatedList->Write();
  myfile->cd();

  TDirectory *sharcdir = myfile->mkdir("SHARC");
  sharcdir->cd();
  sharcList->Write();
  myfile->cd();


  //pCut_array[7]->Write();

  //make a TCanvas of all the individual trific hitmaps:
  auto *c1 = new TCanvas("trif_hitmap_canvas","TRIFIC Hitmaps",800,800);
  c1->Divide(6,4);
  for (int i = 1; i < 25; i++){ //we only do 1-24 because there are only 24 grids in trific
    c1->cd(i);
    trif_hitmap_array[i]->Draw("colz");
  }


  TDirectory *trifdir = myfile->mkdir("TRIFIC");
  trifdir->cd();
  trifList->Write();
  c1->Write();
  sr93PIDCut->Write();
  rb93PIDCut->Write();
  myfile->cd();

  TDirectory *trifHitmapdir = myfile->mkdir("TRIFIC_Hitmaps");
  trifHitmapdir->cd();
  trifHitmapList->Write();
  myfile->cd();

  TDirectory *sharcTig = myfile->mkdir("SHARC-TIGRESS");
  sharcTig->cd();
  sharcTigList->Write();
  myfile->cd();

  TDirectory *trifTig = myfile->mkdir("TRIFIC-TIGRESS");
  trifTig->cd();
  trifTigList->Write();
  myfile->cd();

  TDirectory *sharcTrific = myfile->mkdir("SHARC-TRIFIC");
  sharcTrific->cd();
  sharcTrificList->Write();
  myfile->cd();

  TDirectory *sharcProj = myfile->mkdir("SHARC_Arrays");
  sharcProj->cd();
  sharcArrayList->Write();
  myfile->cd();

  TDirectory *trifBraggGraphs = myfile->mkdir("TRIFIC_Bragg_Graphs");
  trifBraggGraphs->cd();
  trifBraggList->Write();
  myfile->cd();

  TDirectory	*dirKin = myfile->mkdir("Kinematics");
  dirKin->cd(); 
  TGraph *gKin =  sr93_rxn->KinVsTheta();
  gKin->SetName("Sr93_Kinematics");
  gKin->Write();
  myfile->cd();

  /* 
     TDirectory * ggatedir = myfile->mkdir("Gamma_Gated_Focal_Plane");
     ggatedir->cd();
     ggatedList->Write();
     myfile->cd();
  */

  if (bNoTigFilling) std::cout <<"\n\nWarning! Currently skipping all tigress histogram filling!!!!\n\n";fflush(stdout);

  myfile->Write();
  myfile->Close();
}
int main(int argc, char **argv)
{

  SortCode *mysort = new SortCode();

  char const *afile;
  char const *outfile;
  char const *calfile;
  char const *target;
  char const *cutfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0)
  {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();


  // Input-chain-file, output-histogram-file
  if (argc == 1)
  {
    cout << "Insufficient arguments, provide analysis tree" << endl;
    return 0;
  }
  else if (argc == 2)
  {
    afile = argv[1];
    calfile = "CalibrationFile.cal";
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  }
  else if (argc == 3)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = "Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  }
  else if (argc == 4)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  }
  else if (argc == 5)
  {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    target = argv[4];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\nTarget: %s\n", afile, calfile, outfile, target);
    mysort->SortData(afile, calfile, outfile, target);
  }
  else
  {
    printf("Incorrect arguments\n");
    printf("SortData analysis_tree cal_file(optional) out_file(optional)\n");
    return 0;
  }

  

  return 0;
}
