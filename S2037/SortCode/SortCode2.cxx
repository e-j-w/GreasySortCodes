//g++ SortCode.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData
// S1873
// SortCode.cxx
// S. A. Gillespie
// 12/12/2019
// E. J. Williams
// 6/9/2020

#define Sortcode_cxx
#include "SortCode.h"

using namespace std;

Double_t r2d = TMath::RadToDeg();
Double_t d2r = TMath::DegToRad();

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

bool mgate2D(TEmmaHit *hit_one, TCutG *cut)
{
  return cut->IsInside(hit_one->GetPosition().X(), hit_one->GetPosition().Y());
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

double pi = TMath::Pi();
double xp[2] = {-8.0, 8.0};     // PGAC X-position 1D gate, minimum and maximum . Change !!!
double yp[2] = {-15.0, +15.0};  // PGAC T-position 1D gate, minimum and maximum. Change !!!
double tigs3T[2] = {300, 550}; // TIGRESS - S3 Timing. Change !!!
double tigtigT[2] = {-100, 100}; // TIGRESS - TIGRESS Timing. Change !!!

double tempIC;
double tempICArray[4];
double tempIC2D[5][4];
double s3_x_offset = 0.0; // S3 x offset (mm) will need to be recalculated with beam
double s3_y_offset = 0.0; // S3 y offset (mm) will need to be recalculated with beam
double s3_z_offset = 0.0; // S3 z offset (mm) will need to be calculated.
double beta;
double thetalab;
double exc;
double ekin;
double recoiltheta;
double thetacm;
double rekin;

bool suppTig = false;
bool suppAdd = false;
bool s3EnergyDiff = false; //this is for comparing s3 ring vs sector energy differences
bool siicCut = false;      //checking if Si/ic is within the cut we loaded for it

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
  TEmma *emma = 0;
  if (AnalysisTree->FindBranch("TEmma"))
  {
    AnalysisTree->SetBranchAddress("TEmma", &emma);
  }
  else
  {
    cout << "Branch 'TEmma' not found! TEmma variable is NULL pointer" << endl;
  }

  TTigress *tigress = 0;
  if (AnalysisTree->FindBranch("TTigress"))
  {
    AnalysisTree->SetBranchAddress("TTigress", &tigress);
  }
  else
  {
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  // SRIM Files (Need to Change)
  TSRIM *srim_86kr = new TSRIM;
  srim_86kr->ReadEnergyLossFile("kr86_in_sihe.txt");
  TSRIM *srim_94sr = new TSRIM;
  srim_94sr->ReadEnergyLossFile("sr94_in_sihe.txt");
  
  // Set Beam Energy
  double EBeam = 2.6 * 93.915361; //21Na beam
  printf("Beam energy: %f MeV\n", EBeam);

  // Change below with actual target thicknesses, 2.6 micron Si-He
  if (strcmp(target, "one") == 0){
    printf("Using target %s - 2.6 \u03BCm Si-He \n", target); // Assume reaction happens in middle
    EBeam = srim_86kr->GetAdjustedEnergy(EBeam * 1000, 2.6/2., 0.001) / 1000; // GetAdjustedEnergy(energy, target thickness um,step size)
  }else{
    printf("No target specified, not adjusting beam energy.\n");
  }
  printf("Adjusted Beam energy: %f MeV\n", EBeam);

  // Set Reaction
  TReaction *ne21 = new TReaction("sr94", "he4", "n", "zr97", EBeam, 0, true);

  //Defining Pointers
  TEmmaHit *em_hit, *ic_hit, *si_hit, *ssb_hit, *trigger_hit, *anode_hit;
  TTigressHit *tig_hit, *add_hit, *add_hit2;
  //TS3Hit *ring_hit, *sector_hit, *s3hit;
  TVector3 s3pos, recoil_vec;

  //if(s3){
    //s3->SetFrontBackTime(140); // Needed to build S3 pixels properly
    //s3->SetFrontBackTime(1000);
    //s3->SetFrontBackEnergy(0.9);
  //}

  //S3 detector rotation
  //double s3_phi_offset = -90.0*TMath::Pi()/180.; //detector rotated 90deg in chamber
  //s3_phi_offset -= TMath::Pi()/16.; //GRSISort seems to rotate incorrectly by 1 sector width

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
	int 	emma_counter = 0;

	int	counter_ssb = 0;
	int	counter_na22 = 0;

  // recoil beta for doppler correction
  double_t betaDoppler = 0.06562;

  // loop over analysis tree entries.
  printf("\nSorting analysis events...\n");
  for (int jentry = 0; jentry < analentries; jentry++)
  {
    AnalysisTree->GetEntry(jentry);


		if(tigress)
			if(tigress->GetMultiplicity() > 0)
				tig_counter++;-
		if(emma)
			if(emma->GetICMultiplicity() > 0)
				emma_counter++;

		if(tigress && emma){
			if(tigress->GetMultiplicity()>0 && emma->GetICMultiplicity()>0){
				tig_emma_counter++;
			}
		}
    // if tigress event:
    if (tigress)
    {
      // TIGRESS (non addback)
      for (int t = 0; t < tigress->GetMultiplicity(); t++){ // loop over tigress multiplicity (not add-back)
        tig_hit = tigress->GetTigressHit(t);
        if(tig_hit != NULL){
          suppTig = tig_hit->BGOFired();
          if (!suppTig && tig_hit->GetEnergy() > 15){
            tigE->Fill(tig_hit->GetEnergy());
            tigE_ANum->Fill(tig_hit->GetArrayNumber(), tig_hit->GetEnergy());
	    tigE_theta->Fill(tig_hit->GetEnergy(),TMath::RadToDeg()*(tig_hit->GetPosition().Theta()));
          } // suppressor condition
        } // end tig hit condition
      } // end tigress loop

      // TIGRESS addback
      for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){ // loop over add-back events
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
              addE_addE->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy()); // gamma-gamma
              addE_addE->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); // symmetrized
              if (gate1D((add_hit->GetTime() - add_hit2->GetTime()), tigtigT[0], tigtigT[1])){ // tigress timing gate
                addE_addE_tg->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy()); // gamma-gamma
                addE_addE_tg->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
              } // end timing gate
            } // end suppressor condition on second event
          } // end loop over second event 
        } // end suppressor conditon on first event
      } // end loop over first event
    } // end tigress hit condition

    //reset some variables for this new event
    siicCut = false;         //to keep track if the si/ic event is within the cut we made
    //s3EnergyDiff = false;    //this is to compare if the s3 ring vs sector energies are close to each other
    //double ringEnergy = 0.0; //keeping track of total energy in case of split hits
    //double sectorEnergy = 0.0;

    //if (s3)
    //{ // S3 singles
    //  s3->SetMultiHit();
	/*
      for (int i = 0; i < s3->GetRingMultiplicity(); i++)
      {
        ring_hit = s3->GetRingHit(i);
        rings->Fill(ring_hit->GetRing(), ring_hit->GetEnergy());
        ringEnergy += ring_hit->GetEnergy();
        for (int j = 0; j < s3->GetSectorMultiplicity(); j++)
        { //get the ring and sector energies, plot in 2DHisto
          if (s3->GetRingMultiplicity() != 1 || s3->GetSectorMultiplicity() != 1)
            continue; //only look at multiplicity 1 rings and sector hits
          s3_rings_sectors_singles->Fill(ring_hit->GetEnergy(), s3->GetSectorHit(j)->GetEnergy());
          s3_rings_sectors_singlesT->Fill(s3->GetSectorHit(j)->GetTime() - ring_hit->GetTime());
          s3_rings_sectors_singlesTvT->Fill(ring_hit->GetTime()/(1E9*60.), s3->GetSectorHit(j)->GetTime() - ring_hit->GetTime());
          s3_rings_sectors_singlesTvSec->Fill(s3->GetSectorHit(j)->GetSector(), s3->GetSectorHit(j)->GetTime() - ring_hit->GetTime());
          s3_rings_sectors_singlesTvRing->Fill(ring_hit->GetRing(), s3->GetSectorHit(j)->GetTime() - ring_hit->GetTime());
          s3_rings_sectors_singlesTvE->Fill(s3->GetSectorHit(j)->GetEnergy(), s3->GetSectorHit(j)->GetTime() - ring_hit->GetTime());
        }
      }

      for (int j = 0; j < s3->GetSectorMultiplicity(); j++)
      {
        sector_hit = s3->GetSectorHit(j);
        sectors->Fill(sector_hit->GetSector(), sector_hit->GetEnergy());
        sectorEnergy += s3->GetSectorHit(j)->GetEnergy();
      }
      s3_rings_sectors->Fill(ringEnergy, sectorEnergy); //fill this with the sum of rings and sector energies
      if (std::abs(ringEnergy - sectorEnergy) < 100.0)
      { //check if difference in ring and sector energies is less than 100 keV for this event
        s3EnergyDiff = true;
        s3_rings_sectors_gated->Fill(ringEnergy, sectorEnergy);
      }

      
      for (int i = 0; i < s3->GetPixelMultiplicity(); i++){
        s3hit = s3->GetPixelHit(i);
        s3_E->Fill(s3hit->GetEnergy());
        s3pos = s3hit->GetPosition(true);
        s3pos.SetX(s3pos.X() + s3_x_offset);
        s3pos.SetY(s3pos.Y() + s3_y_offset);
        s3pos.SetZ(s3pos.Z() + s3_z_offset);
        s3_E_theta->Fill(s3pos.Theta() * r2d, s3hit->GetEnergy());
        hitmap->Fill(s3pos.X(), s3pos.Y());
        thetalab = s3pos.Theta();
        ekin = s3hit->GetEnergy();
        exc = ne21->GetExcEnergy(ekin * 1e-3, thetalab, 2);
        ne21->SetExcEnergy(exc);
        excE->Fill(exc);
        excE_theta->Fill(s3hit->GetTheta() * r2d, exc);

        if(tigress){
          for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
            add_hit = tigress->GetAddbackHit(t);
            suppAdd = add_hit->BGOFired();
            if (!suppAdd && add_hit->GetEnergy() > 15){
              addT_s3T->Fill(add_hit->GetTime() - s3hit->GetTime());
              if (gate1D((add_hit->GetTime() - s3hit->GetTime()), tigs3T[0], tigs3T[1])){
                //              thetalab = s3pos.Theta();
                //              ekin = s3hit->GetEnergy();
                //              exc = ne21->GetExcEnergy(ekin * 1e-3, thetalab, 2);
                //              ne21->SetExcEnergy(exc);
                //              excE->Fill(exc);
                //              excE_theta->Fill(s3hit->GetTheta() * r2d, exc);
                thetacm = ne21->ConvertThetaLabToCm(thetalab, 2);
                rekin = ne21->GetTLabFromThetaCm(TMath::Pi() - thetacm, 3) * 1e3;
                beta = ne21->AnalysisBeta(rekin * 1e-3, 3);
                recoiltheta = ne21->ConvertThetaCmToLab(thetacm, 3);
                recoil_vec.SetMagThetaPhi(1., recoiltheta, s3pos.Phi() - TMath::Pi());
                addDopp->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                addDopp_ANum->Fill(add_hit->GetArrayNumber(), add_hit->GetDoppler(beta, &recoil_vec));
                addDopp_exc->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                addE_s3_E->Fill(add_hit->GetEnergy(), s3hit->GetEnergy());
                addDopp_s3_E->Fill(add_hit->GetDoppler(beta, &recoil_vec), s3hit->GetEnergy());
              }
            }
          }
          
        }
        
      }
      for (int i = 0; i < s3->GetRingMultiplicity(); i++){
        ring_hit = s3->GetRingHit(i);
        if(tigress){
          for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
            add_hit = tigress->GetAddbackHit(t);
            suppAdd = add_hit->BGOFired();
            if (!suppAdd && add_hit->GetEnergy() > 15){
              //if (gate1D((add_hit->GetTime() - ring_hit->GetTime()), tigs3T[0], tigs3T[1])){
                if(ring_hit->GetEnergy() > 25000){
                  addDopp->Fill(add_hit->GetEnergy());
                }
              //}
            }
          }
        }
      }
    }
*/
    if (emma) { // if EMMA event
      	/*if(emma->GetSSBMultiplicity() > 0)
	//if(emma->GetSSBMultiplicity() > 0 || emma->GetICMultiplicity() > 0 || emma->GetSiMultiplicity() > 0)
	std::cout	<< emma->GetSSBMultiplicity() << "\t" << emma->GetICMultiplicity() << "\t" << emma->GetSiMultiplicity() << "\t" << emma->GetMultiplicity() << "\t" << tigress->GetMultiplicity()
				<< std::endl;*/
      // SSB events
      for (int l = 0; l < emma->GetSSBMultiplicity(); l++){ // loop over SSB events
        ssb_hit = emma->GetSSBHit(l);
        ssbE[ssb_hit->GetDetector()]->Fill(ssb_hit->GetEnergy()); // get detector (left or right), fill energy
        ssbET[ssb_hit->GetDetector()]->Fill(ssb_hit->GetTime() / 1e9, ssb_hit->GetEnergy()); // ssb energy vs time

        // EMMA IC events
        for (int j = 0; j < emma->GetICMultiplicity(); j++){ // loop over IC events
          ic_hit = emma->GetICHit(j);
          ssbICtof->Fill(ssb_hit->GetTime() - ic_hit->GetTime()); // time between ssb and iC? needed?
        } // end IC loop

	// Focal plane silicon events 
        for (int k = 0; k < emma->GetSiMultiplicity(); k++){
          si_hit = emma->GetSiHit(k);
          ssbSItof->Fill(ssb_hit->GetTime() - si_hit->GetTime()); // time between ssb and focal plane silicon
        } // End focal plane silicon loop

	// EMMA triggers sent to TIGRESS
        for (int k = 0; k < emma->GetTriggerMultiplicity(); k++) {
          trigger_hit = emma->GetTriggerHit(k); //EMMA trigger (sent to TIGRESS DAQ)
          ssbEMTtof->Fill(ssb_hit->GetTime() - trigger_hit->GetTime()); // time between ssb and trigger
        } // End EMMA trigger loop
	
	// EMMA TDC events
        for (int m = 0; m < emma->GetTdcMultiplicity(); m++) {
          anode_hit = emma->GetTdcHit(m);
          anodessbT[ssb_hit->GetDetector()]->Fill(ssb_hit->GetTime() - anode_hit->GetTime(), anode_hit->GetSegment());
	break;
        } // end TDC loop
      } // end SSB loop

      // If EMMA + TIGRESS event
      if(tigress){

	// TIGRESS Add-back (coincidence with EMMA)
        for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
          add_hit = tigress->GetAddbackHit(t);

	  // EMMA IC events in coincidence with TIGRESS
          for (int j = 0; j < emma->GetICMultiplicity(); j++){
            ic_hit = emma->GetICHit(j);
            tigICtof->Fill(add_hit->GetTime() - ic_hit->GetTime());
          } // End IC-TIG coincidence event loop

	  // EMMA TDC in coincidence with TIGRESS
          for (int m = 0; m < emma->GetTdcMultiplicity(); m++) {
            anode_hit = emma->GetTdcHit(m);
            tigAnodetof->Fill(add_hit->GetTime() - anode_hit->GetTime());
          } // end EMMA TDC loop in coincidence with TIGRESS
        } // End TIGRESS Addback in coincidence with EMMA 
      } // END TIGRESS event condition

      tempIC = 0; // intialise IC arrays
      for (int k = 0; k < 4; k++)
      {
        tempICArray[k] = 0;
        for (int kk = 0; kk < 5; kk++)
          tempIC2D[kk][k] = 0;
      }

      // EMMA IC singles events
      for (int j = 0; j < emma->GetICMultiplicity(); j++)
      { // Get IC Edep
        ic_hit = emma->GetICHit(j);
        tempICArray[ic_hit->GetSegment() - 1] = ic_hit->GetEnergy();
      } // end emma singles IC event loop

      //   if (emma->GetSiMultiplicity() > 0 && emma->GetSiHit(0)->GetEnergy() > 200) { // Si Gated
      //    if (emma->GetICMultiplicity() == 4 && goodIC(tempICArray) ) {                // IC Gated
      //      if (emma->GetICMultiplicity() == 4) {                // IC Gated
      //if (emma->GetICMultiplicity() == 4 && emma->GetSiMultiplicity() > 0 && emma->GetSiHit(0)->GetEnergy() > 1000)
      //{ // Si + ICGated --GH -- was 200

      // loop over EMMA singles events
      for (int e = 0; e < emma->GetMultiplicity(); e++)
      {
        em_hit = emma->GetEmmaHit(e);

	// EMMA singles Triggers
        for (int one = 0; one < emma->GetTriggerMultiplicity(); one++){
          trigger_hit = emma->GetTriggerHit(one);
          //cout << "hit " << one << ", energy: " << trigger_hit->GetEnergy() << endl;
          if (trigger_hit->GetEnergy() > 200)
            continue;
          tDE->Fill(trigger_hit->GetTime() / pow(10, 9), em_hit->GetTime() - trigger_hit->GetTime()); // get time between EMMA and TIGRESS events
        } // End EMMA singles trigger loop
      } // ENd EMMA singles event loop

      // Loop over IC singles events
      for (int m = 0; m < emma->GetICMultiplicity(); m++) {
	anode_hit = emma->GetICHit(m);

	// Loop over ssb events
	for( int j = 0;j < emma->GetSSBMultiplicity();j++){
		ssb_hit = emma->GetSSBHit(j);
		if(ssb_hit->GetDetector() == 0)
			ssb_anode_dt1->Fill(ssb_hit->GetTime() - anode_hit->GetTime());		
		if(ssb_hit->GetDetector() == 1)
			ssb_anode_dt2->Fill(ssb_hit->GetTime() - anode_hit->GetTime());		
	} // end loop over ssb events
       } // end IC loop

	// loop over focal plane silicon singles
        for (int k = 0; k < emma->GetSiMultiplicity(); k++)
        { // get FP silicon energy
          si_hit = emma->GetSiHit(k);
          siE->Fill(si_hit->GetEnergy());
          siET->Fill(si_hit->GetTime() / 1e9, si_hit->GetEnergy());
        } // End loop over FP silicon singles

        double ICtime = -1;
	// Loop over IC singles
        for (int j = 0; j < emma->GetICMultiplicity(); j++)
        {
          ic_hit = emma->GetICHit(j);
          icE[ic_hit->GetSegment()]->Fill(ic_hit->GetEnergy());
          icN->Fill(ic_hit->GetSegment(), ic_hit->GetEnergy());
          tempIC += ic_hit->GetEnergy();
	  if(ic_hit->GetSegment() == 0)
		ICtime = ic_hit->GetTime();
          if (ic_hit->GetEnergy() < 4000)
            tempICArray[ic_hit->GetSegment() - 1] = ic_hit->GetEnergy(); // IC total energy?
        } // End loop over IC singles

	// Loop over EMMA events
	for (int e = 0; e < emma->GetMultiplicity(); e++)
  	{
    		em_hit = emma->GetEmmaHit(e);
        	xPos->Fill(em_hit->GetPosition().X());
        	yPos->Fill(em_hit->GetPosition().Y());
        	pgac->Fill(em_hit->GetPosition().X(), em_hit->GetPosition().Y());
		// Loop over FP silicon events in coincidence with PGAC
		for(int k=0;k<emma->GetSiMultiplicity();k++){
   	       		si_hit = emma->GetSiHit(k);
			if(ne22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)) // PID cut 1
				pgac_ne22->Fill(em_hit->GetPosition().X(), em_hit->GetPosition().Y()); // PID gated 2D PGAC position
			if(na22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)) // PID cut 2
				pgac_na22->Fill(em_hit->GetPosition().X(), em_hit->GetPosition().Y()); // PID gated 2D PGAC position
		} // End loop over FP silicon in coincidence with PGAC
	} // End loop over PGAC events

        if (tempICArray[0] != 0 && tempICArray[1] != 0 && tempICArray[2] != 0 && tempICArray[3] != 0) // if hit in all IC segments
        { // build 2D IC related histograms
          icSum->Fill(tempIC);
          iC0V1->Fill(tempICArray[0], tempICArray[1]);
          iC0V2->Fill(tempICArray[0], tempICArray[2]);
          iC0V3->Fill(tempICArray[0], tempICArray[3]);
          iC1V2->Fill(tempICArray[1], tempICArray[2]);
          iC1V3->Fill(tempICArray[1], tempICArray[3]);
          iC2V3->Fill(tempICArray[2], tempICArray[3]);
	  iC0VSum->Fill(tempICArray[0],tempIC);
	  // If TIGRESS + EMMA event
          if(tigress){
		// Loop over TIGRESS Addback events in coincidence with EMMA
          	for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
			if(add_hit->GetEnergy() < 50)
				continue;
            		add_hit = tigress->GetAddbackHit(t);
			if(add_hit->GetTime() - ICtime < 800 || add_hit->GetTime() - ICtime > 1050)
				continue;
			if(add_hit->GetDoppler(betaDoppler) > 1250 && add_hit->GetDoppler(betaDoppler) < 1300)
				iC0VSum_1274->Fill(tempICArray[0],tempIC);
			if(ne22_cut_ic->IsInside(tempICArray[0],tempIC))
				addDopp_ne22_ic->Fill(add_hit->GetDoppler(beta));
			if(na22_cut_ic->IsInside(tempICArray[0],tempIC))
				addDopp_na22_ic->Fill(add_hit->GetDoppler(beta));
		} // End loop over TIGRESS addback events in coincidence with EMMA
	  } // End Tigress event condition

	  // Loop over FP silicon events
          for (int k = 0; k < emma->GetSiMultiplicity(); k++)
          { // build hists for IC vs Silicon
            si_hit = emma->GetSiHit(k);
            icSumVSi->Fill(si_hit->GetEnergy(), tempIC);
            icSumVicPlusSi->Fill(tempIC+si_hit->GetEnergy(),tempIC);
            ic0VSi->Fill(si_hit->GetEnergy(), tempICArray[0]);
            ic1VSi->Fill(si_hit->GetEnergy(), tempICArray[1]);
            ic2VSi->Fill(si_hit->GetEnergy(), tempICArray[2]);
            ic3VSi->Fill(si_hit->GetEnergy(), tempICArray[3]);

            /*if(s3){
              if(s3->GetPixelMultiplicity() > 0){
                icSumVicPlusSi_s3gate->Fill(tempIC+si_hit->GetEnergy(),tempIC);
              }
            }*/
				
	     if(ne22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC))
		iC0VSum_ne22->Fill(tempICArray[0],tempIC);
	     if(na22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC))
		iC0VSum_na22->Fill(tempICArray[0],tempIC);

	    // If TIGRESS in coincidence with EMMA
            if(tigress){
		// Loop over TIGRESS addback events
              for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++){
                add_hit = tigress->GetAddbackHit(t);
		if(add_hit->GetDoppler(betaDoppler) > 1250 && add_hit->GetDoppler(betaDoppler) < 1300)
			ICvsSiICSum_1274->Fill(tempIC+si_hit->GetEnergy(),tempIC);
                suppAdd = add_hit->BGOFired();
                if (!suppAdd && add_hit->GetEnergy() > 15){
                  //printf("energy: %f\n",add_hit->GetDoppler(betaDoppler));
                  if(add_hit->GetDoppler(betaDoppler) > 1254){
                    if(add_hit->GetDoppler(betaDoppler) < 1294){
                      icSumVicPlusSi1274Gate->Fill(tempIC+si_hit->GetEnergy(),tempIC);
                    } // End <1294 keV gate
                  } // End > 1254 keV gate
                } // End suppressor condition
              } // End Add back event loop
            } // End Tigress event condition
          } // End FP silicon event loop
        } // End all IC hit condition

	for( int s=0;s<emma->GetSSBMultiplicity();s++){
		counter_ssb++;
	}
	for( int e=0;e<emma->GetSiMultiplicity();e++){
	        si_hit = emma->GetSiHit(e);
		tempIC = 0;
		for(int j=0;j<emma->GetICMultiplicity();j++){
			ic_hit = emma->GetICHit(j);
			tempIC += ic_hit->GetEnergy();
		}
		emma_dE_E->Fill(si_hit->GetEnergy(),tempIC);
		if(na22_cut->IsInside(tempIC+si_hit->GetEnergy(),tempIC))
			counter_na22++;				
		if ( tigress ) {
			//std::cout	<< "BOO " << tigress->GetMultiplicity() << std::endl;
			for( int t =0; t< tigress->GetAddbackMultiplicity(); t++){
				add_hit	= tigress->GetAddbackHit(t);
				suppAdd = add_hit->BGOFired();
				addemmatof->Fill(add_hit->GetTime() - si_hit->GetTime());
       				addE_tof->Fill((add_hit->GetTime() - si_hit->GetTime()), add_hit->GetEnergy());
            			addemmatofdet->Fill(add_hit->GetArrayNumber(),add_hit->GetTime() - si_hit->GetTime());
				emma_tig_dt_E->Fill(si_hit->GetTime()-add_hit->GetTime(),add_hit->GetDoppler(betaDoppler));
				if(add_hit->GetTime() - si_hit->GetTime() > 820 && add_hit->GetTime() - si_hit->GetTime() < 1000){
					double beta = 0.0859;
					if(!suppAdd && add_hit->GetEnergy()>50){
						addDopp_emma->Fill(add_hit->GetDoppler(beta));
                				if(add_hit->GetArrayNumber()>48){
                  					addDopp_emma_lampshade->Fill(add_hit->GetDoppler(beta));
                				}else{
                 					addDopp_emma_corona->Fill(add_hit->GetDoppler(beta));
                				}
						if(ne22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)){
							addDopp_ne22->Fill(add_hit->GetDoppler(beta));
							add_ne22vt->Fill(add_hit->GetEnergy(),add_hit->GetPosition().Theta()*TMath::RadToDeg());
							for(int tt=t+1;tt<tigress->GetAddbackMultiplicity();tt++){
										add_hit2 = tigress->GetAddbackHit(tt);
										suppAdd = add_hit2->BGOFired();
										if(!suppAdd && add_hit2->GetEnergy()>50){
											if(add_hit->GetTime() - add_hit2->GetTime() > -100 && add_hit->GetTime() - add_hit2->GetTime() < 100){
												addDopp_ne22_gg->Fill(add_hit->GetDoppler(beta),add_hit2->GetDoppler(beta));
												addDopp_ne22_gg->Fill(add_hit2->GetDoppler(beta),add_hit->GetDoppler(beta));
											}
										}
									}
								}
								if(na22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)){
									addDopp_na22->Fill(add_hit->GetDoppler(beta));
									add_na22vt->Fill(add_hit->GetEnergy(),add_hit->GetPosition().Theta()*TMath::RadToDeg());
									for(int tt=t+1;tt<tigress->GetAddbackMultiplicity();tt++){
										add_hit2 = tigress->GetAddbackHit(tt);
										suppAdd = add_hit2->BGOFired();
										if(!suppAdd && add_hit2->GetEnergy()>50){
											if(add_hit->GetTime() - add_hit2->GetTime() > -100 && add_hit->GetTime() - add_hit2->GetTime() < 100){
												addDopp_na22_gg->Fill(add_hit->GetDoppler(beta),add_hit2->GetDoppler(beta));
												addDopp_na22_gg->Fill(add_hit2->GetDoppler(beta),add_hit->GetDoppler(beta));
											}
										}
									}
								}
                if (mg22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)){
                  addDopp_mg22->Fill(add_hit->GetDoppler(beta));
                  if(add_hit->GetArrayNumber()>48){
                    addDopp_mg22_lampshade->Fill(add_hit->GetDoppler(beta));
                  }else{
                    addDopp_mg22_corona->Fill(add_hit->GetDoppler(beta));
                  }
                  for(int tt=t+1;tt<tigress->GetAddbackMultiplicity();tt++){
										add_hit2 = tigress->GetAddbackHit(tt);
										suppAdd = add_hit2->BGOFired();
										if(!suppAdd && add_hit2->GetEnergy()>50){
											if(add_hit->GetTime() - add_hit2->GetTime() > -100 && add_hit->GetTime() - add_hit2->GetTime() < 100){
												addDopp_mg22_gg->Fill(add_hit->GetDoppler(beta),add_hit2->GetDoppler(beta));
												addDopp_mg22_gg->Fill(add_hit2->GetDoppler(beta),add_hit->GetDoppler(beta));
											}
										}
									}            
                }
							}
              //unsuppressed
              if(add_hit->GetEnergy()>50){
                if (mg22_cut->IsInside(tempIC + si_hit->GetEnergy(),tempIC)){
                  if(add_hit->GetArrayNumber()>48){
                    addDopp_mg22_lampshade_unsupp->Fill(add_hit->GetDoppler(beta));
                  }else{
                    addDopp_mg22_corona_unsupp->Fill(add_hit->GetDoppler(beta));
                  }
                }
              }
              

						}

					}

				}
			}

			if( s3 && emma->GetICMultiplicity()>0){


					double ICtime = -1;
 		      for (int j = 0; j < emma->GetMultiplicity(); j++)
   		   	{
	    	    em_hit = emma->GetEmmaHit(j);
            for (int i = 0; i < s3->GetPixelMultiplicity(); i++)
            {
              s3hit = s3->GetPixelHit(i);
						  s3pos = s3hit->GetPosition(true);
              s3pos.SetX(s3pos.X() + s3_x_offset);
              s3pos.SetY(s3pos.Y() + s3_y_offset);
              s3pos.SetZ(s3pos.Z() + s3_z_offset);
              s3emmatof->Fill(s3hit->GetTime() - em_hit->GetTime());

						if(s3hit->GetTime() - em_hit->GetTime() > 450 && s3hit->GetTime() - em_hit->GetTime() < 620){

              thetalab = s3pos.Theta();
              ekin = s3hit->GetEnergy();
              exc = ne21->GetExcEnergy(ekin * 1e-3, thetalab, 2); //Energy conversion from keV to MeV

							emma_s3_exc->Fill(exc);
							emma_s3_exc_ring->Fill(s3hit->GetRing(),exc);
			
							if(tigress){
								for(int t=0;t<tigress->GetAddbackMultiplicity();t++){
									add_hit = tigress->GetAddbackHit(t);
									suppAdd = add_hit->BGOFired();
									if(!suppAdd && add_hit->GetEnergy() > 50 && add_hit->GetTime() - em_hit->GetTime() > 850 && add_hit->GetTime() - em_hit->GetTime() < 1020){
										emma_s3_exc_addback->Fill(add_hit->GetDoppler(beta),exc);
									}
								}
							}
						}

					}
					}
		
	//	}

	}

         /* if (s3 && tigress)
          { // EMMA + S3 + TIGRESS coincidences

            for (int i = 0; i < s3->GetPixelMultiplicity(); i++)
            {
              s3hit = s3->GetPixelHit(i);
              for (int t = 0; t < tigress->GetAddbackMultiplicity(); t++)
              {
                add_hit = tigress->GetAddbackHit(t);
                suppAdd = add_hit->BGOFired();
                if (!suppAdd && add_hit->GetEnergy() > 15)
                {
                  if (gate1D((add_hit->GetTime() - s3hit->GetTime()), tigs3T[0], tigs3T[1]))
                  {
                    s3pos = s3hit->GetPosition(true);
                    s3pos.SetX(s3pos.X() + s3_x_offset);
                    s3pos.SetY(s3pos.Y() + s3_y_offset);
                    s3pos.SetZ(s3pos.Z() + s3_z_offset);

                    thetalab = s3pos.Theta();
                    ekin = s3hit->GetEnergy();
                    exc = ne21->GetExcEnergy(ekin * 1e-3, thetalab, 2); //Energy conversion from keV to MeV
                    ne21->SetExcEnergy(exc);
                    thetacm = ne21->ConvertThetaLabToCm(thetalab, 2);
                    rekin = ne21->GetTLabFromThetaCm(TMath::Pi() - thetacm, 3) * 1e3;
                    beta = ne21->AnalysisBeta(rekin * 1e-3, 3); //get laboratory beta of the recoil
                    recoiltheta = ne21->ConvertThetaCmToLab(thetacm, 3);
                    recoil_vec.SetMagThetaPhi(1., recoiltheta, s3pos.Phi() - TMath::Pi());

                    s3emmatof->Fill(s3hit->GetTime() - em_hit->GetTime());

                    if (tofGate((s3hit->GetTime() - em_hit->GetTime())))
                    { // TOF gate and
                      excE_tg->Fill(exc);
                      excE_theta_tg->Fill(s3hit->GetTheta() * r2d, exc);
                      addDopp_tg->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                      addDopp_ANum_tg->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
                      addDopp_exc_tg->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                      s3_E_gated->Fill(s3hit->GetEnergy());
                      s3_E_theta_gated->Fill(s3pos.Theta() * r2d, s3hit->GetEnergy());
                      hitmap_time_gated->Fill(s3pos.X(), s3pos.Y());

                      //TIGRESS-TIGRESS-S3 time gated
                      for (int t2 = t+1; t2 < tigress->GetAddbackMultiplicity(); t2++)
                      {
                        add_hit2 = tigress->GetAddbackHit(t2);
                        suppAdd = add_hit2->BGOFired();
                        if (gate1D((add_hit->GetTime() - add_hit2->GetTime()), tigtigT[0], tigtigT[1]))
                          {
                            addE_addE_tofg->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy());
                            addE_addE_tofg->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); //symmetrized
                            addDopp_addDopp_tg->Fill(add_hit->GetDoppler(beta, &recoil_vec), add_hit2->GetDoppler(beta, &recoil_vec));
                            addDopp_addDopp_tg->Fill(add_hit2->GetDoppler(beta, &recoil_vec), add_hit->GetDoppler(beta, &recoil_vec)); //symmetrized
                          }
                      }*/

                      /*if (emma->GetICMultiplicity() == 4 && emma->GetSiMultiplicity() > 0 && emma->GetSiHit(0)->GetEnergy() > 1000 && add_hit->GetDoppler(beta, &recoil_vec) > 340 && add_hit->GetDoppler(beta, &recoil_vec) < 350)
                      { //add in the Si vs IC with the 350 keV gamma gate
                        for (int k = 0; k < emma->GetSiMultiplicity(); k++)
                        {

                          icSumVSi_gated_350->Fill(emma->GetSiHit(k)->GetEnergy(), tempIC);
                        }
                      }

                      
                      if (gate1D(em_hit->GetPosition().X(), xp[0], xp[1]) && gate1D(em_hit->GetPosition().Y(), yp[0], yp[1]))
                      { // PGAC position gated
                        excE_PIDG->Fill(exc);
                        excE_theta_PIDG->Fill(s3hit->GetTheta() * r2d, exc);
                        addDopp_PIDG->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                        addDopp_ANum_PIDG->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
                        addDopp_exc_PIDG->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                      }
                      */
                    }
                    /*else if (BGtofGate((s3hit->GetTime() - em_hit->GetTime())))
                    { // TOF gate
                      excE_bg->Fill(exc);
                      excE_theta_bg->Fill(s3hit->GetTheta() * r2d, exc);
                      addDopp_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                      addDopp_ANum_bg->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
                      addDopp_exc_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                      if (CutG_loaded)
                      {
                        //if(mgate2D(em_hit,massGate)) {
                        if (siicCut && s3EnergyDiff)
                        {
                          excE_PIDG_bg->Fill(exc);
                          excE_theta_PIDG_bg->Fill(s3hit->GetTheta() * r2d, exc);
                          addDopp_PIDG_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                          addDopp_ANum_PIDG_bg->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
                          addDopp_exc_PIDG_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                        }
                      }
                      else
                      {
                        if (gate1D(em_hit->GetPosition().X(), xp[0], xp[1]) && gate1D(em_hit->GetPosition().Y(), yp[0], yp[1]))
                        { // PGAC position gated + TOF gate
                          excE_PIDG_bg->Fill(exc);
                          excE_theta_PIDG_bg->Fill(s3hit->GetTheta() * r2d, exc);
                          addDopp_PIDG_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec));
                          addDopp_ANum_PIDG_bg->Fill(add_hit->GetArrayNumber(), add_hit->GetEnergy());
                          addDopp_exc_PIDG_bg->Fill(add_hit->GetDoppler(beta, &recoil_vec), exc);
                        }
                      }
                    }*/
           /*       }
                }
              }
            }
          } // if S3 & tigress
        }  // emma
      //}     //si
    }  */     //if emma
    if (jentry % 10000 == 0)
      cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
  } // analysis tree
  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;

	std::cout	<< analentries << " events, " << tig_emma_counter << " containing TIGRESS + EMMA, " << tig_counter << " containing TIGRESS, " << emma_counter << " containing EMMA"
						<< std::endl;
	std::cout	<< "Na22: " << counter_na22 << "\t SSB: " << counter_ssb 
			<< std::endl;

  cout << "Writing histograms to " << outfile << endl;

  TFile *myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory *tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  myfile->cd();

  TDirectory *s3dir = myfile->mkdir("S3");
  s3dir->cd();
  s3List->Write();
  myfile->cd();

  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();

  TDirectory *tigs3dir = myfile->mkdir("TIGRESS-S3");
  tigs3dir->cd();
  tigs3List->Write();
  myfile->cd();

  TDirectory *emmadir = myfile->mkdir("EMMA");
  emmadir->cd();
  emmaList->Write();
  TDirectory *emmaicdir = emmadir->mkdir("IC");
  emmaicdir->cd();
  icList->Write();
  emmadir->cd();
  myfile->cd();

  TDirectory *s3emmadir = myfile->mkdir("EMMA-S3");
  s3emmadir->cd();
  emmas3List->Write();
  myfile->cd();

  TDirectory *tigemmadir = myfile->mkdir("TIGRESS-EMMA");
  tigemmadir->cd();
  tigemmaList->Write();
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

  TDirectory *s3pidgatedir = myfile->mkdir("PID_Gated_S3");
  s3pidgatedir->cd();
  PIDgatedListS3->Write();
  myfile->cd();

  /* 
        TDirectory * ggatedir = myfile->mkdir("Gamma_Gated_Focal_Plane");
        ggatedir->cd();
        ggatedList->Write();
        myfile->cd();*/

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
