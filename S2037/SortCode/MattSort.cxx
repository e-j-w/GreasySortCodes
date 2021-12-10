// g++ MattSort.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SortData_Matt
// S1873
// SortCode_Matt.cxx
// M. Williams
// 17/12/2019

#define MattSort_cxx
#include "MattSort.h"


using namespace std;

Double_t r2d = TMath::RadToDeg();

bool gate1D(Double_t value, Double_t min, Double_t max)
{
  if (min < value && value < max)
    return true;
  else
    return false;
}

double tigtigT[2] = {-100, 100};

double pi = TMath::Pi();
double Beta;
double thetalab; 
double excite; 
double ekin;
double recoiltheta;
double thetacm;
double rekin;
double betaDoppler = 0.07005;

bool suppAdd = false;

void MattSort::SortData(char const * afile, char const * calfile, char const * outfile, char const * target = "NULL", char const * cutfile = "NULL") {
  
	// rootfile branch structure
	typedef struct {double xpos, ypos, multi;} PGAC;
	PGAC pgac;
	typedef struct {double energy, time, multi;} SILICON;
	SILICON silicon;
	typedef struct {double eseg[5], esum, multi;} IONCHAMBER;
	IONCHAMBER ic;
	typedef struct {double energy[2], time[2], multi;} SSB;
	SSB sb;
	typedef struct {double time, multi;} TOF;
	TOF tof;
	typedef struct {double time;} TRIG;
	TRIG trig;
	typedef struct {double rawE, dopE, time, det, core, beta;} GAMMA;
	GAMMA gamma;
	typedef struct {double rawE, dopE, time, det, core, beta;} GAMMA2;
	GAMMA2 gamma2;

	// define output rootfile and branches
	TFile * myfile = new TFile(outfile, "RECREATE");
	// SSB
	TTree * ssb = new TTree("ssb","SSB Singles");
	ssb->Branch("sb",&sb,"energy[2]/D:time[2]/D:multi/D");
	// EMMA singles
  	TTree * emma = new TTree("emma","EMMA Singles");
  	emma->Branch("pgac",&pgac,"xpos/D:ypos/D:multi/D");
	emma->Branch("silicon",&silicon,"energy/D:time/D:multi/D");
	emma->Branch("ic",&ic,"eseg[5]/D:esum/D:multi/D");
	//emma->Branch("tof",&tof,"time/D:multi/D");
	emma->Branch("trig",&trig,"time/D");
	// TIGRESS Singles
	TTree * tig = new TTree("tig","TIGRESS Singles");
  	tig->Branch("gamma",&gamma,"rawE/D:dopE/D:time/D:det/D:core/D:beta/D");
	tig->Branch("gamma2",&gamma2,"rawE/D:dopE/D:time/D:det/D:core/D:beta/D");
	// EMMA-TIGRESS coincidences
	TTree * tigemma = new TTree("tigemma","TIGRESS-EMMA Coincidences");
	tigemma->Branch("pgac",&pgac,"xpos/D:ypos/D:multi/D");
	tigemma->Branch("silicon",&silicon,"energy/D:time/D:multi/D");
	tigemma->Branch("ic",&ic,"eseg[5]/D:esum/D:multi/D");
	tigemma->Branch("tof",&tof,"time/D:multi/D");
	tigemma->Branch("gamma",&gamma,"rawE/D:dopE/D:time/D:det/D:core/D:beta/D");
	tigemma->Branch("gamma2",&gamma2,"rawE/D:dopE/D:time/D:det/D:core/D:beta/D");
	
  // Open Analysis Tree
  TFile * analysisfile = new TFile(afile, "READ");   
  if (!analysisfile->IsOpen()) {
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }
 
  printf("File %s opened\n", afile);
  TChain * AnalysisTree = (TChain * ) analysisfile->Get("AnalysisTree");
  long analentries = AnalysisTree->GetEntries();
  const char * testval = "NULL";

  // Checks for branches and sets pointers
  // Set EMMA hits
  TEmma * emma_data = 0;
  if (AnalysisTree->FindBranch("TEmma")) {
    AnalysisTree->SetBranchAddress("TEmma", & emma_data);
  } else {
    cout << "Branch 'TEmma' not found! TEmma variable is NULL pointer" << endl;
  }
  TTigress * tig_data = 0;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tig_data);
  } else {
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }
  // Get SRIM files
  // SRIM Files (Need to Change)
  TSRIM *srim_86kr = new TSRIM;
  srim_86kr->ReadEnergyLossFile("kr86_in_sihe.txt");
  TSRIM *srim_94sr = new TSRIM;
  srim_94sr->ReadEnergyLossFile("sr94_in_sihe.txt");

  double EBeam = 2.6 * 93.915361; // 94Sr beam

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
  TReaction *ne21 = new TReaction("sr94", "he4", "n", "zr97", EBeam, 0, true); // Not needed?

  // Define Hit Pointers
  TEmmaHit * em_hit, * si_hit, * ic_hit, * ssb_hit, * trigger_hit;
  TTigressHit * tig_hit, * add_hit, * add_hit2;
  TVector3 pos, recoil_vec;

  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile); // Read calibration file
  printf("\nSorting analysis events...\n");

  for (int jentry = 0; jentry < analentries; jentry++) { // loop over events in analysis tree
  
	AnalysisTree->GetEntry(jentry);

	if (tig_data)	{
		for (int t = 0; t < tig_data->GetAddbackMultiplicity(); t++) { // Loop over tigress events
        		tig_hit = tig_data->GetAddbackHit(t);
        		//tig_hit = tig_data->GetTigressHit(t);
       			// if(tig_hit != NULL){
          		suppAdd = tig_hit->BGOFired();
          		if (!suppAdd && tig_hit->GetEnergy() > 15) { // suppressor condition?
				gamma.rawE = tig_hit->GetEnergy(); // get doppler corrected gamma energy
				gamma.dopE = tig_hit->GetDoppler(betaDoppler); // get doppler corrected gamma energy
				gamma.time = tig_hit->GetTime();
				gamma.det = tig_hit->GetDetector();
				tig->Fill(); // Fill TIGRESS singles tree
			} // End Suppressor Conditon
    // }
		} // End Loop over TIGRESS singles Events
	} // End tig hit

	if (emma_data) { // EMMA

		// SSB data
		sb.energy[2] = {0};
		sb.time[2] = {0};
		sb.multi = emma_data->GetSSBMultiplicity(); // SSB Multi
		for (int l = 0; l < emma_data->GetSSBMultiplicity(); l++) { // Get SSB hits
        		ssb_hit = emma_data->GetSSBHit(l);
			sb.energy[ssb_hit->GetDetector()] = ssb_hit->GetEnergy(); // ssb time (function 'GetDetector' is indexed starting from 1)
        		sb.time[ssb_hit->GetDetector()] = ssb_hit->GetTime()/1e9; // ssb time
			ssb->Fill(); // Fill SSB tree
      		}
		
		// EMMA event multiplicities
		pgac.multi = emma_data->GetMultiplicity(); // Get EMMA Multiplicity
		silicon.multi = emma_data->GetSiMultiplicity(); // Get EMMA Si Multi
		ic.multi = emma_data->GetICMultiplicity(); // Get IC Multi
		tof.multi = emma_data->GetTriggerMultiplicity(); // Get Trigger Multi

		// PGAC hits
		for (int e = 0; e < emma_data->GetMultiplicity(); e++) { // loop over emma events
			em_hit = emma_data->GetEmmaHit(e);
			pgac.xpos = em_hit->GetPosition().X(); // Get PGAC xpos
      			pgac.ypos = em_hit->GetPosition().Y(); // Get PGAC ypos
     			// Trigger
			for (int m = 0; m < emma_data->GetTriggerMultiplicity(); m++) {
				trig.time = 0.;
        			trigger_hit = emma_data->GetTriggerHit(m);
        			//if (trigger_hit->GetEnergy() > 200) continue;
           				//tDE->Fill(trigger_hit->GetTime() / pow(10, 9), em_hit->GetTime() - trigger_hit->GetTime()); // get time between EMMA and TIGRESS events
				//tof.time = em_hit->GetTime() - trigger_hit->GetTime();
				trig.time = trigger_hit->GetTime();
   			}

		}
		// Silicon data
		silicon.energy = 0;
		silicon.time = 0;
		for (int k = 0; k < emma_data->GetSiMultiplicity(); k++) {
			si_hit = emma_data->GetSiHit(k);
			silicon.energy = si_hit->GetEnergy(); // get FP silicon energy
            		silicon.time = si_hit->GetTime() / 1e9; // get FP silicon time		
		}

		// IC data
		ic.eseg[5] = {0};
		ic.esum = 0;
		for (int j = 0; j < emma_data->GetICMultiplicity(); j++) { // IC multi = how many segments have a hit
            		ic_hit = emma_data->GetICHit(j);
            		ic.eseg[ic_hit->GetSegment()]=ic_hit->GetEnergy(); // get ic segement energy
            		ic.esum += ic_hit->GetEnergy(); // get ic sum energy
          	}
		emma->Fill(); // Fill EMMA singles tree		
	}

	if (tig_data && emma_data) { // TIGRESS & S3 hit condition
	// EMMA event multiplicities
		pgac.multi = emma_data->GetMultiplicity(); // Get EMMA Multiplicity
		silicon.multi = emma_data->GetSiMultiplicity(); // Get EMMA Si Multi
		ic.multi = emma_data->GetICMultiplicity(); // Get IC Multi
		tof.multi = emma_data->GetTriggerMultiplicity(); // Get Trigger Multi
		// PGAC hits
		for (int e = 0; e < emma_data->GetMultiplicity(); e++) { // loop over emma events
			em_hit = emma_data->GetEmmaHit(e);
			pgac.xpos = em_hit->GetPosition().X(); // Get PGAC xpos
      			pgac.ypos = em_hit->GetPosition().Y(); // Get PGAC ypos
			// Trigger
			for (int m = 0; m < emma_data->GetTriggerMultiplicity(); m++) {
        			trigger_hit = emma_data->GetTriggerHit(m);
        			//if (trigger_hit->GetEnergy() > 200) continue;
        			//tDE->Fill(trigger_hit->GetTime() / pow(10, 9), em_hit->GetTime() - trigger_hit->GetTime()); // get time between EMMA and TIGRESS events
			  	
      			}
		// IC data
		ic.eseg[5] = {0};
		ic.esum = 0;
		for (int j = 0; j < emma_data->GetICMultiplicity(); j++) { // IC multi = how many segments have a hit
        		ic_hit = emma_data->GetICHit(j);
        		ic.eseg[ic_hit->GetSegment()]=ic_hit->GetEnergy(); // get ic segement energy
        		ic.esum += ic_hit->GetEnergy(); // get ic sum energy
      		}
		// TIGRESS data
		for (int t = 0; t < tig_data->GetAddbackMultiplicity(); t++) { // Loop over tigress events
       			add_hit = tig_data->GetAddbackHit(t);
          		suppAdd = add_hit->BGOFired();
          		if (!suppAdd && add_hit->GetEnergy() > 15) { // suppressor condition?
				gamma.dopE = add_hit->GetDoppler(betaDoppler); // get doppler corrected gamma energy
				gamma.time = add_hit->GetTime();
				gamma.rawE = add_hit->GetEnergy();
				gamma.det = add_hit->GetDetector();
				tof.time = add_hit->GetTime() - em_hit->GetTime(); // TIGRESS - EMMA Time
				/*// gamma-gamma coinc
					for (int t2 = t+1; t2 < tig_data->GetAddbackMultiplicity(); t2++) {
            					add_hit2 = tig_data->GetAddbackHit(t2);
              					suppAdd = add_hit2->BGOFired();
              					if (gate1D((add_hit->GetTime() - add_hit2->GetTime()), tigtigT[0], tigtigT[1])) {
                					gamma2.dopE = add_hit2->GetDoppler(Beta, & recoil_vec);
						    	gamma2.time = add_hit2->GetTime();
							gamma2.rawE = add_hit2->GetEnergy();
							tigemma->Fill(); // Fill S3-TIGRESS-EMMA coincidence tree
              					} // end gamma-gamma time gate
            				}*/ // end 2nd gamma loop
				tigemma->Fill(); // Fill TIGRESS-EMMA coincidence tree 
			} // End Suppressor Conditon
		} // End Loop over TIGRESS Events
		} // End loop over emma events
	} // End TIGRESS & EMMA hit condition
	
  	if (jentry % 10000 == 0)
      		cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush; // Event counter
  } // End Analysis Tree Loop

  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;

  myfile->cd();
  myfile->Write();
  myfile->Close();
}

int main(int argc, char ** argv) {


  MattSort * mysort = new MattSort();

  char const * afile;
  char const * outfile;
  char const * calfile;
  char const * target;
  char const * cutfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if (grsi_path.length() > 0) {
    grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);
  TParserLibrary::Get()->Load();

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
    cout << "Insufficient arguments, provide analysis tree" << endl;
    return 0;
  } else if (argc == 2) {
    afile = argv[1];
    //calfile = "CalibrationFile.cal";
    //calfile = "CalibrationFile2020Nov20.cal";
    calfile = "CalibrationFile_Dec19.cal";
    outfile = "HistFiles/Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  } else if (argc == 3) {
    afile = argv[1];
    calfile = argv[2];
    outfile = "HistFiles/Histograms.root";
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  } else if (argc == 4) {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    printf("Analysis file: %s#define Sortcode_cxx\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
    mysort->SortData(afile, calfile, outfile);
  } else if (argc == 5) {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    target = argv[4];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\nTarget: %s\n", afile, calfile, outfile, target);
    mysort->SortData(afile, calfile, outfile, target);
  } else if (argc == 6) {
    afile = argv[1];
    calfile = argv[2];
    outfile = argv[3];
    target = argv[4];
    cutfile = argv[5];
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\nTarget: %s\nCuts File: %s\n", afile, calfile, outfile, target, cutfile);
    mysort->SortData(afile, calfile, outfile, target, cutfile);
  } else if (argc > 6) {
    printf("Doh! Too many arguments\n");
    return 0;
  }

  return 0;
}
