//g++ SortCode.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -lTTigress -lTEmma -lTGRSIDetector -lTGenericDetector -o SortData
// S1641 
// SortCode.cxx
// S. A. Gillespie
// 02/09/2019

#define Sortcode_cxx

#include "SortCode.h"

using namespace std;

Double_t r2d = TMath::RadToDeg();

bool tofGate(Double_t Time){
	if (Time < 950 && Time > 500) return true;
	else return false;
}

//ToF background gate
bool BGtofGate(Double_t Time){
	if (Time < 0 && Time > -1500) return true;
	if (Time < 2500 && Time > 1000) return true;
	return false;
}

bool mgate1D(Double_t Position, Double_t min, Double_t max) {
  	if (min < Position && Position < max) return true;
	else return false;
}


bool mgate2D(TEmmaHit* hit_one, TCutG* cut){
   	return cut->IsInside(hit_one->GetPosition().X(), hit_one->GetPosition().Y());
}

bool ggate1D(Double_t Energy, Double_t min, Double_t max) {
  	if (min < Energy && Energy < max) return true;
	else return false;
}

bool loadCutG(char const *cutfile) {
	TFile *cuts = new TFile(cutfile,"READ");
	massGate[0] = (TCutG*)cuts->Get("mass_1");
	massGate[1] = (TCutG*)cuts->Get("mass_2");
	//icGate[0] = (TCutG*)cuts->Get("ic_1");
	//icGate[1] = (TCutG*)cuts->Get("ic_2");
	return true;
}

double xMin[5] = {-100,-100,-100,-100,-100};
double xMax[5] = {100,100,100,100,100};
double gMin[5] = {788,659,718,0,0};
double gMax[5] = {800,666,723,819,4096};
double tempIC;
double tempICArray[4];
double tempIC2D[5][4];
double beta;
bool CutG_loaded = false;
bool suppTig = false;
bool suppAdd = false;
void SortCode::SortData(char const * afile, char const * calfile, char const * outfile, char const * target = "NULL", char const * cutfile = "NULL") {
  Initialise();
  TFile * analysisfile = new TFile(afile, "READ"); //Opens Analysis Trees
  if (!analysisfile->IsOpen()) {
    printf("Opening file %s failed, aborting\n", afile);
    return;
  }
  printf("File %s opened\n", afile);
  TChain * AnalysisTree = (TChain * ) analysisfile->Get("AnalysisTree");
  long analentries = AnalysisTree->GetEntries();
  const char * testval = "NULL";
  if (strcmp(cutfile, testval) != 0) {
    printf("Cuts Loaded from %s \n",cutfile);
    loadCutG(cutfile);
    CutG_loaded = true;
  } else printf("Couldn't find 2D cuts... 1D gates will be used instead.\n");

  // Checks for branches and sets pointers 
  TEmma * emma = 0;
  if (AnalysisTree->FindBranch("TEmma")) {
    AnalysisTree->SetBranchAddress("TEmma", & emma);
  } else {
    cout << "Branch 'TEmma' not found! TEmma variable is NULL pointer" << endl;
  }

  TTigress * tigress = 0;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
  } else {
    cout << "Branch 'TTigress' not found! TTigress variable is NULL pointer" << endl;
  }

  if (strcmp(target, "thick") == 0) {
    printf("Using %s target - 1300 ug/cm^3 \n",target);
    beta = 0.07005; // Target ? 2.4 MeV/u
//   beta = 0.0693; // Target 3 2.7 MeV/u
//   beta = 0.0664; //Target 1
  } else if (strcmp(target, "middle") == 0) {
    printf("Using %s target - 550 ug/cm^3 \n",target);
    beta = 0.07005;
  } else if (strcmp(target, "thin") == 0) {
    printf("Using %s target - 250 ug/cm^3 \n",target);
    beta = 0.07005;
  } else {
    printf("No Target Specified Assuming Thick Target.\n");
    beta = 0.07005;
  }

  //Defining Pointers
  TEmmaHit * em_hit, * ic_hit, * si_hit, *ssb_hit, *trigger_hit;
  TTigressHit * tig_hit, *add_hit, *add_hit2;
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);

  printf("\nSorting analysis events...\n");
	for (int jentry = 0; jentry < analentries; jentry++) { // Loop over analysis tree entries

		AnalysisTree->GetEntry(jentry);
		if (tigress) {	// If TIGRESS events
			for(int t=0; t<tigress->GetMultiplicity(); t++){ // Loop over TIGRESS singles 
				tig_hit = tigress->GetTigressHit(t);
				suppTig = tig_hit->BGOFired();
				if(!suppTig && tig_hit->GetEnergy()>15) {
					tigE->Fill(tig_hit->GetEnergy());
					tigDop->Fill(tig_hit->GetDoppler(beta));
					// begin CRN
					//tigE_tigDop->Fill(tig_hit->GetEnergy(),tig_hit->GetDoppler(beta));
					//end CRN
					tigE_ANum->Fill(tig_hit->GetArrayNumber(),tig_hit->GetEnergy());
					tigDop_ANum->Fill(tig_hit->GetArrayNumber(),tig_hit->GetDoppler(beta));
					tigE_theta->Fill(tig_hit->GetPosition().Theta()*r2d,tig_hit->GetEnergy());
					tigDop_theta->Fill(tig_hit->GetPosition().Theta()*r2d,tig_hit->GetDoppler(beta));
					for(int bgoInd=0; bgoInd < tigress->GetBGOMultiplicity(); bgoInd++){
            if((tigress->GetBGO(bgoInd).GetDetector() == tig_hit->GetDetector()) && (tigress->GetBGO(bgoInd).GetEnergy() > 0.)){
              tigT_bgoT_supp->Fill(tig_hit->GetCfd() - tigress->GetBGO(bgoInd).GetCfd());
            }
          }
				} // END suppressor conditon
				if(tig_hit->GetEnergy() > 15){

          bgo_mult->Fill(tigress->GetBGOMultiplicity());
          for(int bgoInd=0; bgoInd < tigress->GetBGOMultiplicity(); bgoInd++){
            bgo_det->Fill(tigress->GetBGO(bgoInd).GetDetector());
            if((tigress->GetBGO(bgoInd).GetDetector() == tig_hit->GetDetector()) && (tigress->GetBGO(bgoInd).GetEnergy() > 0.)){
              tigT_bgoT->Fill(tig_hit->GetCfd() - tigress->GetBGO(bgoInd).GetCfd());
            }
          }
        }
			} // End Tigress event loop
			

			for(int t=0; t<tigress->GetAddbackMultiplicity(); t++){ 
				add_hit = tigress->GetAddbackHit(t);
				suppAdd = add_hit->BGOFired();
				if(!suppAdd && add_hit->GetEnergy()>15) {
					addE->Fill(add_hit->GetEnergy());
					addDop->Fill(add_hit->GetDoppler(beta));
					addE_ANum->Fill(add_hit->GetArrayNumber(),add_hit->GetEnergy());
					addDop_ANum->Fill(add_hit->GetArrayNumber(),add_hit->GetDoppler(beta));
					addE_theta->Fill(add_hit->GetPosition().Theta()*r2d,add_hit->GetEnergy());
					addDop_theta->Fill(add_hit->GetPosition().Theta()*r2d,add_hit->GetDoppler(beta));
					//TIGRESS-TIGRESS addback
					for (int t2 = t+1; t2 < tigress->GetAddbackMultiplicity(); t2++){
						add_hit2 = tigress->GetAddbackHit(t2);
						suppAdd = add_hit2->BGOFired();
						if (!suppAdd && add_hit2->GetEnergy() > 15){
							addT_addT->Fill(add_hit->GetTime() - add_hit2->GetTime());
							addE_addE->Fill(add_hit->GetEnergy(),add_hit2->GetEnergy()); // gamma-gamma
							addE_addE->Fill(add_hit2->GetEnergy(),add_hit->GetEnergy()); // symmetrized
						} // end suppressor condition on second event
					} // end loop over second event 
				} // End suppressor condition
			} // End Addback event loop
		} // End TIGRESS conditoon

		// Temp IC array:
		tempIC=0;
		for(int k = 0; k<4;k++){
				tempICArray[k]=0;
				for(int kk = 0; kk<5;kk++)tempIC2D[kk][k]=0;
		}
   		
   		/*for( int e=0;e<emma->GetSiMultiplicity();e++){
			si_hit = emma->GetSiHit(e);	
			if ( tigress ) {
				//std::cout	<< "BOO " << tigress->GetMultiplicity() << std::endl;
				for( int t =0; t< tigress->GetAddbackMultiplicity(); t++){
					add_hit	= tigress->GetAddbackHit(t);
					suppAdd = add_hit->BGOFired();
					if(!suppTig && tig_hit->GetEnergy()>15){ // Suppressor and Array hit condition
						addemmatof->Fill(add_hit->GetTime() - si_hit->GetTime());
		       				addE_tof->Fill((add_hit->GetTime() - si_hit->GetTime()), add_hit->GetEnergy());
		       				addDop_tof->Fill((add_hit->GetTime() - si_hit->GetTime()),add_hit->GetDoppler(beta));
					}		
				}
			}
		}*/
   		
		if(emma){

			for(int l = 0; l < emma->GetSSBMultiplicity(); l++) { // Loop over SSB singles
					ssb_hit = emma->GetSSBHit(l);
					ssbE[ssb_hit->GetDetector()]->Fill(ssb_hit->GetEnergy()*0.001);
					ssbET[ssb_hit->GetDetector()]->Fill(ssb_hit->GetTime()/1e9,ssb_hit->GetEnergy()*0.001);
				} // End SSB loop

			for (int e = 0; e < emma->GetMultiplicity(); e++) { // Loop over EMMA events

				

				for (int k = 0; k < emma->GetSiMultiplicity(); k++) { // Focal plane silicon singles
					si_hit = emma->GetSiHit(k);
					siE->Fill(si_hit->GetEnergy());
					siET->Fill(si_hit->GetTime()/1e9,si_hit->GetEnergy());
				} // End FP silicon event loop

				em_hit = emma->GetEmmaHit(e);
				xPos->Fill(em_hit->GetPosition().X());
				yPos->Fill(em_hit->GetPosition().Y());
				xPosT->Fill(em_hit->GetTime()/1e9,em_hit->GetPosition().X());
				yPosT->Fill(em_hit->GetTime()/1e9,em_hit->GetPosition().Y());
				pgac->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());				

				for(int t=0; t<tigress->GetMultiplicity(); t++){ // Loop over TIGRESS singles 
					tig_hit = tigress->GetTigressHit(t);
					suppTig = tig_hit->BGOFired();
					if(!suppTig && tig_hit->GetEnergy()>15) {
						// begin CRN
						//tigE_tigDop->Fill(tig_hit->GetEnergy(),tig_hit->GetDoppler(beta));
						//end CRN
						
					} // END suppressor conditon
				} // End Tigress event loop
		

				for (int one = 0; one < emma->GetTriggerMultiplicity(); one++) {
					trigger_hit = emma->GetTriggerHit(one);
					if(trigger_hit->GetEnergy() > 200)continue;
					tDE->Fill(trigger_hit->GetTime()/pow(10,9),em_hit->GetTime()-trigger_hit->GetTime());
				} // End Loop over EMMA triggers

				// Loop over IC events
				for (int j = 0; j < emma->GetICMultiplicity(); j++) { 
					ic_hit = emma->GetICHit(j);
					icE[ic_hit->GetSegment()]->Fill(ic_hit->GetEnergy());
					icN->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
					tempIC += ic_hit->GetEnergy();
					if(ic_hit->GetEnergy()<4000) tempICArray[ic_hit->GetSegment()-1] = ic_hit->GetEnergy();
				} // End loop over IC events

				if(tempICArray[0]!=0 && tempICArray[1]!=0 && tempICArray[2]!=0 && tempICArray[3]!=0){
					icSum->Fill(tempIC);
					iC0V1->Fill(tempICArray[0],tempICArray[1]);
					iC0V2->Fill(tempICArray[0],tempICArray[2]);
					iC0V3->Fill(tempICArray[0],tempICArray[3]);
					iC1V2->Fill(tempICArray[1],tempICArray[2]);
					iC1V3->Fill(tempICArray[1],tempICArray[3]);
					iC2V3->Fill(tempICArray[2],tempICArray[3]);
					// Loop over Si events
					for(int k=0; k < emma->GetSiMultiplicity(); k++) {
						si_hit = emma->GetSiHit(k);
						icSumVSi->Fill(si_hit->GetEnergy(),tempIC);
						ic0VSi->Fill(si_hit->GetEnergy(),tempICArray[0]);
						ic1VSi->Fill(si_hit->GetEnergy(),tempICArray[1]);
						ic2VSi->Fill(si_hit->GetEnergy(),tempICArray[2]);
						ic3VSi->Fill(si_hit->GetEnergy(),tempICArray[3]);
					} // End loop over Si events
				} // End IC hit condition

				// If TIGRESS + EMMA  
				if (tigress){
					for(int t=0; t<tigress->GetMultiplicity(); t++){  // Loop over TIGRESS events
						tig_hit = tigress->GetTigressHit(t);
						suppTig = tig_hit->BGOFired();
						if(!suppTig && tig_hit->GetEnergy()>15){ // Suppressor and Array hit condition
							tigDop_ANum_if_emma->Fill(tig_hit->GetArrayNumber(),tig_hit->GetDoppler(beta));
							tigemmatof->Fill(tig_hit->GetTime() - em_hit->GetTime());
							tigE_tof->Fill((tig_hit->GetTime() - em_hit->GetTime()),tig_hit->GetEnergy()); 
							tigDop_tof->Fill((tig_hit->GetTime() - em_hit->GetTime()),tig_hit->GetDoppler(beta));
							if( tofGate((tig_hit->GetTime() - em_hit->GetTime())) ){ // EMMA-TIGRESS timing gate condition
								tigE_tg->Fill(tig_hit->GetEnergy());
								if(tig_hit->GetArrayNumber()<48){
									addDop_tg->Fill(tig_hit->GetDoppler(beta));
								} 
								tigE_xPos_tg->Fill(em_hit->GetPosition().X(),tig_hit->GetEnergy()); 
								tigDop_xPos_tg->Fill(em_hit->GetPosition().X(),tig_hit->GetDoppler(beta)); 
								pgac_tg->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
								if(CutG_loaded){
									for(int i=0; i < 2; i++) {
										if( mgate2D( em_hit, massGate[i])) {
												tigE_massG[i]->Fill(tig_hit->GetEnergy());
												tigDop_massG[i]->Fill(tig_hit->GetDoppler(beta));
										} // end mass gate condition
									} // end loop
								}else{ // if cut not loaded
									for(int i=0; i < 2; i++) {
										if( mgate1D( em_hit->GetPosition().X(), xMin[i], xMax[i])) { // mass gate xpos condition
												tigE_massG[i]->Fill(tig_hit->GetEnergy());
												tigDop_massG[i]->Fill(tig_hit->GetDoppler(beta));
										} // end mass gate condition
									} // end loop 
								} // End If cut loaded
								for(int i=0; i < 4; i++) {
									if( ggate1D(tig_hit->GetDoppler(beta), gMin[i], gMax[i]) ) { // gamma gate condition
										pgac_gGate[i]->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
									} // end gamma gate condition
								} // end loop
							}else if( BGtofGate((tig_hit->GetTime() - em_hit->GetTime())) ){ // if does not satisfy suppressor condition? 
								if(tig_hit->GetArrayNumber()<48){
									addDop_tg_bg->Fill(tig_hit->GetDoppler(beta));
								}
								tigE_xPos_bg->Fill(em_hit->GetPosition().X(),tig_hit->GetEnergy()); 
								tigDop_xPos_bg->Fill(em_hit->GetPosition().X(),tig_hit->GetDoppler(beta)); 
								pgac_bg->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
								if(CutG_loaded) { // If graphical cut loaded
									for(int i=0; i < 2; i++) {
										if( mgate2D( em_hit, massGate[i]) ) { // if mass gate satisfied
											tigE_massG_bg[i]->Fill(tig_hit->GetEnergy());
											tigDop_massG_bg[i]->Fill(tig_hit->GetDoppler(beta));
										} // end mass gate condition
									} // end loop
								}else{
									for(int i=0; i < 2; i++) {
										if( mgate1D( em_hit->GetPosition().X(), xMin[i], xMax[i])) {
											tigE_massG_bg[i]->Fill(tig_hit->GetEnergy());
											tigDop_massG_bg[i]->Fill(tig_hit->GetDoppler(beta));
										} // End Mass gate condition
									} // End loop
								} // end else
								for(int i=0; i < 4; i++) {
									if( ggate1D(tig_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
										pgac_gGate_bg[i]->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
									} // End gamma codition
								} // End loop
							} // End BGO condition
						} // End Tigress loop
					} // End Tigress non add back event loop

					// Add back event loop (coincidence with EMMA)
					for(int t=0; t<tigress->GetAddbackMultiplicity(); t++){
						add_hit = tigress->GetAddbackHit(t);
						suppAdd = add_hit->BGOFired();
						if(!suppAdd && add_hit->GetEnergy()>15){
							addemmatof->Fill(add_hit->GetTime() - em_hit->GetTime());
							addE_tof->Fill((add_hit->GetTime() - em_hit->GetTime()),add_hit->GetEnergy()); 
							addDop_tof->Fill((add_hit->GetTime() - em_hit->GetTime()),add_hit->GetDoppler(beta));
							for(int i=0; i < 2; i++) {
								if(CutG_loaded){
									if( mgate2D( em_hit, massGate[i])) {
										addDop_mtof[i]->Fill(add_hit->GetTime() - em_hit->GetTime());
									} // end 2D mass gate
								}
							} // End loop
							if( tofGate((add_hit->GetTime() - em_hit->GetTime())) ){
								addE_xPos_tg->Fill(em_hit->GetPosition().X(),add_hit->GetEnergy()); 
								addDop_xPos_tg->Fill(em_hit->GetPosition().X(),add_hit->GetDoppler(beta)); 
								if(CutG_loaded){
									for(int i=0; i < 2; i++) {
										if( mgate2D( em_hit, massGate[i])) {
											addE_massG[i]->Fill(add_hit->GetEnergy());
											addDop_massG[i]->Fill(add_hit->GetDoppler(beta));
										} // End mass gate condition
									} // End loop	
								}else{
									for(int i=0; i < 2; i++) {
										if( mgate1D( em_hit->GetPosition().X(), xMin[i], xMax[i])) {
											addE_massG[i]->Fill(add_hit->GetEnergy());
											addDop_massG[i]->Fill(add_hit->GetDoppler(beta));
										} // End emma mass gate condition
									} // End loop
								}
								for(int i=0; i < 4; i++) {
									if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
										pgac_aGate[i]->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
									} // End Gamma gate condition 
								} // End loop
								if(CutG_loaded) {
									if( mgate2D( em_hit, massGate[0])) {
										for (int j = 0; j < emma->GetICMultiplicity(); j++) { 
											ic_hit = emma->GetICHit(j);
												icN_tg->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
											for(int i=0; i < 4; i++) {
												if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
													icN_aGate_tg[i]->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
												tempIC2D[i][ic_hit->GetSegment()-1] = ic_hit->GetEnergy();
												}
													if(tempIC2D[i][0]!=0 && tempIC2D[i][1]!=0 && tempIC2D[i][2]!=0 && tempIC2D[i][3]!=0){
												iC0V1_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][1]);
												iC0V2_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][2]);
												iC0V3_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][3]);
												iC1V2_aGate_tg[i]->Fill(tempIC2D[i][1],tempIC2D[i][2]);
												iC1V3_aGate_tg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
												iC2V3_aGate_tg[i]->Fill(tempIC2D[i][2],tempIC2D[i][3]);
											for(int k=0; k < emma->GetSiMultiplicity(); k++) {
												si_hit = emma->GetSiHit(k);
												icSumVSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),(tempIC2D[i][0]+tempIC2D[i][1]+tempIC2D[i][2]+tempIC2D[i][3]));
												ic0VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][0]);
												ic1VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][1]);
												ic2VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][2]);
												ic3VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][3]);
											} // end FP silicon event loop
													} // End IC hit condition
												} // End loop
										} // End loop over IC events
									} // End mass gate condition
									else if( mgate1D( em_hit->GetPosition().X(), xMin[0], xMax[0])) {
										for (int j = 0; j < emma->GetICMultiplicity(); j++) { 
											ic_hit = emma->GetICHit(j);
												icN_tg->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
											for(int i=0; i < 4; i++) {
												if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
													icN_aGate_tg[i]->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
												tempIC2D[i][ic_hit->GetSegment()-1] = ic_hit->GetEnergy();
												}
													if(tempIC2D[i][0]!=0 && tempIC2D[i][1]!=0 && tempIC2D[i][2]!=0 && tempIC2D[i][3]!=0){
												iC0V1_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][1]);
												iC0V2_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][2]);
												iC0V3_aGate_tg[i]->Fill(tempIC2D[i][0],tempIC2D[i][3]);
												iC1V2_aGate_tg[i]->Fill(tempIC2D[i][1],tempIC2D[i][2]);
												iC1V3_aGate_tg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
												iC2V3_aGate_tg[i]->Fill(tempIC2D[i][2],tempIC2D[i][3]);
											for(int k=0; k < emma->GetSiMultiplicity(); k++) {
												si_hit = emma->GetSiHit(k);
												icSumVSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),(tempIC2D[i][0]+tempIC2D[i][1]+tempIC2D[i][2]+tempIC2D[i][3]));
												ic0VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][0]);
												ic1VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][1]);
												ic2VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][2]);
												ic3VSi_aGate_tg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][3]);
											} // End loop over silicon events
													} // End IC hit condition
												} // End loop
										} // End loop over IC events
									} // End X-pos gate condition 	 
								} // End Cut load condition
							}else if( BGtofGate((add_hit->GetTime() - em_hit->GetTime())) ){
									addE_xPos_bg->Fill(em_hit->GetPosition().X(),add_hit->GetEnergy()); 
									addDop_xPos_bg->Fill(em_hit->GetPosition().X(),add_hit->GetDoppler(beta)); 
									if(CutG_loaded) {
										for(int i=0; i < 2; i++) {
												if( mgate2D( em_hit, massGate[i])) {
														addE_massG_bg[i]->Fill(add_hit->GetEnergy());
														addDop_massG_bg[i]->Fill(add_hit->GetDoppler(beta));
												} // End 2D mass gate
										} // End loop
									} else {
										for(int i=0; i < 2; i++) {
												if( mgate1D( em_hit->GetPosition().X(), xMin[i], xMax[i])) {
														addE_massG_bg[i]->Fill(add_hit->GetEnergy());
														addDop_massG_bg[i]->Fill(add_hit->GetDoppler(beta));
												} // End 1D mass gate
										} // End loop
									} // end else
									for(int i=0; i < 4; i++) {
										if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
												pgac_aGate_bg[i]->Fill(em_hit->GetPosition().X(),em_hit->GetPosition().Y());
										} // End gamma gate
									} // end loop

									if(CutG_loaded) {
									if( mgate2D( em_hit, massGate[0])) {
										for (int j = 0; j < emma->GetICMultiplicity(); j++) {
											ic_hit = emma->GetICHit(j);
											icN_bg->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
											for(int i=0; i < 4; i++) {
												if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
													icN_aGate_bg[i]->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
													tempIC2D[i][ic_hit->GetSegment()-1] = ic_hit->GetEnergy();
												}
												if(tempIC2D[i][0]!=0 && tempIC2D[i][1]!=0 && tempIC2D[i][2]!=0 && tempIC2D[i][3]!=0){
													iC0V1_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][1]);
													iC0V2_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][2]);
													iC0V3_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][3]);
													iC1V2_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][2]);
													iC1V3_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
													iC2V3_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
													for(int k=0; k < emma->GetSiMultiplicity(); k++) {
														si_hit = emma->GetSiHit(k);
														icSumVSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),(tempIC2D[i][0]+tempIC2D[i][1]+tempIC2D[i][2]+tempIC2D[i][3]));
														ic0VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][0]);
														ic1VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][1]);
														ic2VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][2]);
														ic3VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][3]);
													} // End loop over Si events
												} // End IC hit Condtion
											} // End loop
										} // End IC event loop
									} // End 2D mass gate
										} // End load Gcut	
								else if( mgate1D( em_hit->GetPosition().X(), xMin[0], xMax[0])) {
									for (int j = 0; j < emma->GetICMultiplicity(); j++) {
										ic_hit = emma->GetICHit(j);
										icN_bg->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
										for(int i=0; i < 4; i++) {
											if( ggate1D(add_hit->GetDoppler(beta), gMin[i], gMax[i]) ) {
												icN_aGate_bg[i]->Fill(ic_hit->GetSegment(),ic_hit->GetEnergy());
												tempIC2D[i][ic_hit->GetSegment()-1] = ic_hit->GetEnergy();
											}
											if(tempIC2D[i][0]!=0 && tempIC2D[i][1]!=0 && tempIC2D[i][2]!=0 && tempIC2D[i][3]!=0){
												iC0V1_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][1]);
												iC0V2_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][2]);
												iC0V3_aGate_bg[i]->Fill(tempIC2D[i][0],tempIC2D[i][3]);
												iC1V2_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][2]);
												iC1V3_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
												iC2V3_aGate_bg[i]->Fill(tempIC2D[i][1],tempIC2D[i][3]);
												for(int k=0; k < emma->GetSiMultiplicity(); k++) {
													si_hit = emma->GetSiHit(k);
													icSumVSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),(tempIC2D[i][0]+tempIC2D[i][1]+tempIC2D[i][2]+tempIC2D[i][3]));
													ic0VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][0]);
													ic1VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][1]);
													ic2VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][2]);
													ic3VSi_aGate_bg[i]->Fill(si_hit->GetEnergy(),tempIC2D[i][3]);
												} // End silicon event loop
											} // End IC hit condition
										} // End loop
									} // End loop over IC events
								} // End X-pos gate condition
							}
						} // End suppressor condition
					} // End addback events loop
				} // End tigress condition
			} // End EMMA event loop
		} // End EMMA hit condition

    if (jentry % 10000 == 0){
			cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << analentries << ", " << 100 * jentry / analentries << "% complete" << "\r" << flush;
		}

  } // End analysis event loop
  cout << "Entry " << analentries << " of " << analentries << ", 100% complete" << endl;
  cout << "Event sorting complete" << endl;

  cout << "Writing histograms to " << outfile << endl;
  
  addDop_tg_bgs->Add(addDop_tg,addDop_tg_bg,1.0,-0.15);

  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  TDirectory * emmadir = myfile->mkdir("EMMA");
  emmadir->cd();
  emmaList->Write();
  TDirectory * emmaicdir = emmadir->mkdir("IC");
  emmaicdir->cd();
  icList->Write();
  emmadir->cd();
  myfile->cd();
  TDirectory * tigdir = myfile->mkdir("TIGRESS");
  tigdir->cd();
  tigList->Write();
  addList->Write();
  myfile->cd();
  TDirectory *tigtigdir = myfile->mkdir("TIGRESS-TIGRESS");
  tigtigdir->cd();
  tigtigList->Write();
  myfile->cd();
  TDirectory *tigbgodir = myfile->mkdir("TIGRESS-BGO");
  tigbgodir->cd();
  tigbgoList->Write();
  myfile->cd();
  TDirectory * tigemmadir = myfile->mkdir("TIGRESS-EMMA");
  tigemmadir->cd();
  tigemmaList->Write();
  addemmaList->Write();
  myfile->cd();
  TDirectory * tgatedir = myfile->mkdir("ToF_Gated_TIGRESS");
  tgatedir->cd();
  timegatedList->Write();
  myfile->cd();
  TDirectory * mgatedir = myfile->mkdir("Mass_Gated_TIGRESS");
  mgatedir->cd();
  massgatedList->Write();
  myfile->cd();
  TDirectory * ggatedir = myfile->mkdir("Gamma_Gated_Focal_Plane");
  ggatedir->cd();
  ggatedList->Write();
  TDirectory * icggatedir = ggatedir->mkdir("IC");
  icggatedir->cd();
  ggatedicList->Write(); 
  ggatedir->cd();
  TDirectory * icsiggatedir = ggatedir->mkdir("IC-Si");
  icsiggatedir->cd();
  ggatedICSiList->Write();
  ggatedir->cd();
  myfile->cd();
  myfile->Write();
  myfile->Close();
}

int main(int argc, char ** argv) {

  SortCode * mysort = new SortCode();

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
    calfile = "CalibrationFile.cal";
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
    printf("Analysis file: %s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);
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
