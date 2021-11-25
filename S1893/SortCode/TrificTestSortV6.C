{
    string FileName="TrifOut.root";
    double factor=1;
    
    
    // Loading analysis tree etc (this is all done for you in grsiproof)
 ///////// FILL IN A FEW MANUAL ANALYSIS TREE FILE NAMES ///////// 
    TChain* DataChain = new TChain("AnalysisTree","AnalysisTree");

//     TChannel::ReadCalFile("CalFiles/Historic.cal");
//     DataChain->Add("Data/OldRun/analysis51989_000.root");FileName="OldRun.root";factor=1;
//     DataChain->Add("Data/OldRun/analysis51994_000.root");FileName="OldRun94.root";factor=0.1;
    
    
//     TChannel::ReadCalFile("CalFiles/CalibrationFile.cal");
//     DataChain->Add("Data/analysis55528_000.root");FileName="Run28.root";factor=0.1;
//     DataChain->Add("Data/analysis55529_000.root");FileName="Run29.root";factor=1;
//     DataChain->Add("Data/analysis55530_000.root");FileName="Run30.root";factor=1;
//     DataChain->Add("Data/analysis55535_000.root");FileName="Run35.root";factor=1;
    
    TChannel::ReadCalFile("/tig/pterodon_data3/S1893/NewTrifCalFile.cal");
    DataChain->Add("/tig/pterodon_data3/S1893/AnalysisTrees/analysis55540_000.root");FileName="Run40.root";factor=1;
//     DataChain->Add("Data/analysis55542_000.root");FileName="Run42.root";factor=1;
    DataChain->Add("/tig/pterodon_data3/S1893/AnalysisTrees//analysis55543_000.root");FileName="Run43.root";factor=1;
//     DataChain->Add("Data/analysis55544_000.root");FileName="Run44.root";factor=0.1;
    FileName="Run4X.root";
    
//     TChannel::ReadCalFile("CalFiles/NewTrificCal.cal");
	/// AAA ~2.5 epA x100 attenuator
//     DataChain->Add("Data/BadDAQ/anallysis55502_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55503_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55504_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55505_000.root");
	
	/// BBB ~10 epA epA x100 attenuator
//     DataChain->Add("Data/BadDAQ/anallysis55507_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55508_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55509_000.root");
// 	
// 	// CCC >20 epA x10 attenuator 
//     DataChain->Add("Data/BadDAQ/anallysis55510_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55511_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55512_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55513_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55514_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55518_000.root");
//     DataChain->Add("Data/BadDAQ/anallysis55519_000.root");

    // Set address of branch
    TTrific *trific=0;
    if(DataChain->FindBranch("TTrific"))DataChain->SetBranchAddress("TTrific",&trific);else return;
	cout<<endl<<"LOADED BRANCH"<<endl;
	
    // Load things and set things for class statics
    TTrific::SetSharc();
//     TVector3 gridnormal=TTrific::fNormalGridVec;
    double fAngle = (60./180.)*TMath::Pi();
    TVector3 gridnormal=TVector3(0,-TMath::Cos(fAngle),TMath::Sin(fAngle));
    
    double fInitialSpacingCart=TTrific::fInitialSpacingCart;
    double fSpacingCart=TTrific::fSpacingCart;

    int FirstGridN=1;


	TFile out(FileName.c_str(),"RECREATE");
    out.cd();
		TH2F* RawSingles=new TH2F("RawSingles","RawSingles;Grid Number;Charge",24,0,24,1000,0,3000);
		TH2F* RawSinglesLow=new TH2F("RawSinglesLow","RawSinglesLow;Grid Number;Charge",24,0,24,1000,0,200);
		TH1D* dTimeSingles=new TH1D("dTimeSingles","Time Difference Single Grids;dT (ns);cunts",200,-1000,1000);
		TH1D* dTlong=new TH1D("dTlong","Time Difference Single Grids;dT (ns);cunts",1000,0,10000);
		TH2F* SingleMult=new TH2F("SingleMult","Single Grid Multiplicity;Grid Number;Multiplicity",24,0,24,10,0,10);
		
		TH1D* dTlong5er=new TH1D("dTlong5er","Time Difference Single Grids>5;dT (ns);cunts",1000,0,10000);
		
		out.mkdir("kValues");
		out.cd("kValues");
			TH2F* SingleskValue=new TH2F("SingleskValue","Single Grid kValue;Grid Number;log10(kValue)",32,0,32,20,-1,4);
			TH2F* XkValue=new TH2F("XkValue","X Grid kValue;Segment Number;log10(kValue)",12,0,12,20,-1,4);
			TH2F* YkValue=new TH2F("YkValue","Y Grid kValue;Segment Number;log10(kValue)",12,0,12,20,-1,4);
		out.cd();
		
		out.mkdir("TimeGated");
		out.cd("TimeGated");	
			TH1D* dTGated=new TH1D("dTGated","Time Difference Time Gated Single Grids;dT (ns);cunts",200,-1000,1000);
			TH2F* GatedSingles=new TH2F("GatedSingles","GatedSingles;Grid Number;Charge",24,0,24,1000,0,3000);
		out.cd();
		
		out.mkdir("GrigGated");
		out.cd("GrigGated");	

			TH2F* TrigCount=new TH2F("TrigCount","TrigCount;Grid Number of Trigger;Grid Number",24,0,24,24,0,24);
			TH2F* Singles10=new TH2F("Singles10","Singles10;Grid Number;Charge",24,0,24,1000,0,3000);
			TH2F* Singles12=new TH2F("Singles12","Singles12;Grid Number;Charge",24,0,24,1000,0,3000);
			
			
		out.cd();
        
        out.mkdir("Correlated");
		out.cd("Correlated");	
            TH2F* ContinuitySingles=new TH2F("ContinuitySingles","ContinuitySingles;Grid Number;Charge",24,0,24,1000,0,3000);

			TH2F* ID_RangeZ=new TH2F("ID_RangeZ","ID_RangeZ",24,0,24,1000,0,1000);
			TH2F* ID_Range4=new TH2F("ID_Range4","ID_Range4",24,0,24,1000,0,1000);
			TH2F* ID_Range2=new TH2F("ID_Range2","ID_Range2",24,0,24,1000,0,1000);
			TH2F* ID_Z2=new TH2F("ID_Z2","ID_Z2",1000,0,1000,1000,0,1000);
			TH2F* ID_Z4=new TH2F("ID_Z4","ID_Z4",1000,0,1000,1000,0,1000);
			
			
		out.cd();
        
		TH2F* XSingles=new TH2F("XSingles","XSingles",12,0,12,1000,0,3000);
		TH2F* YSingles=new TH2F("YSingles","YSingles",12,0,12,1000,0,3000);
		TH2F* XYMult=new TH2F("XYMult","XYMult",12,0,12,12,0,12);
		TH2F* XYMAPRAW=new TH2F("XYMAPRAW","XYMAPRAW",12,0,12,12,0,12);
		TH2F* XYmap_single=new TH2F("XYmap_single","XYmap_single",12,0,12,12,0,12);
		
        
        TH1D* dTimeY6=new TH1D("dTimeY5","dTimeY5;dT (ns);cunts",200,-1000,1000);
        TH1D* YN=new TH1D("YN","YN",13,0,13);
        TH1D* YN6=new TH1D("YN6","YN6",13,0,13);
        
        TH2D* EY=new TH2D("EY","EY",13,0,13,2000,0,2000);
        
        TH2D* EX=new TH2D("EX","EX",13,0,13,2000,0,2000);
        
        
        TH2D* dTEY6=new TH2D("dTEY6","dTEY6;dT (ns);Energy",200,-1000,1000,2000,0,2000);
        TH2D* dTEX6=new TH2D("dTEX6","dTEX6;dT (ns);Energy",200,-1000,1000,2000,0,2000);
        
        
        TH2D* dT6=new TH2D("dT6","dT6;dT (ns);Grid",200,-1000,1000,24,0,24);
        
    gROOT->cd();
	
    vector<int> Rfact={0,1};
    for(int i=2;i<24;i++){
        Rfact.push_back(Rfact[i-1]*i);
        if(i==3)Rfact[i]=Rfact[i-1];
        if(i==5)Rfact[i]=Rfact[i-1];
    }
    Rfact[3]=0;
    Rfact[5]=0;
    
    Int_t nentries = DataChain->GetEntries();
	
	vector< vector<TTrificHit*>> TimeCorrelatedSingles;
	vector< double > dE2;
	vector< double > dE4;
	vector< double > dEMax;
	vector< int > Range;
	vector< int > Contin;
	
    // Loop over events
    for(int jentry = 0; jentry<nentries*factor ; jentry++){
		
		// Loop frontmatter
		TimeCorrelatedSingles.clear();
		dE2.clear();
		dE4.clear();
		dEMax.clear();
		Range.clear();
		Contin.clear();
        
        
        DataChain->GetEntry(jentry);
        if(jentry%1000 == 0){
            cout << setiosflags(ios::fixed) << std::setprecision(2) << 100.0*jentry/nentries << " % complete.\r" << flush;
        }
		
		// holder for raw multiplicity
		int Ncount[24] = { 0 };
		
        Double_t T6=0;
        
        double tE12=0;
        double tE10=0;
        
		// loop over trific single grids
        for(int h=0;h<trific->GetSingMultiplicity();h++){
			// Fetch  hit
            TTrificHit* hit=trific->GetTrificSingHit(h);
			int d=hit->GetDetector();
            double tE=hit->GetEnergy();
			
            if(d==6)T6=hit->GetTime();
            
			// Fill multiplicity counter
			if(d<24&&tE>0)Ncount[d]++;
			
			//fill raw energy
            RawSingles->Fill(d,tE);
			RawSinglesLow->Fill(d,tE);
			// Fill singles kValue graph
			SingleskValue->Fill(hit->GetDetector(),log10(hit->GetKValue()));
			
            // Skip the noise
            if(tE<10)continue;
            
            if(d==12)tE12=tE;
            if(d==10)tE10=tE;
            
            if(d==3){cout<<"endl grid3 error"<<endl;continue;}
            if(d==5){cout<<"endl grid5 error"<<endl;continue;}
            
            
			// Loop over all other single grids to compare time
			for(int k=h+1;k<trific->GetSingMultiplicity();k++){
				TTrificHit* hitb=trific->GetTrificSingHit(k);
				dTimeSingles->Fill(hitb->GetTime()-hit->GetTime());
				dTlong->Fill(hitb->GetTime()-hit->GetTime());
				if(d>5&&hitb->GetDetector()>5)dTlong5er->Fill(hitb->GetTime()-hit->GetTime());
			}
			
			// Build TimeCorrelatedSingles sub-events
			bool UnUsed=true;
			for(unsigned int k=0;k<TimeCorrelatedSingles.size();k++){ //loop over existing sub-events
				for(unsigned int j=0;j<TimeCorrelatedSingles[k].size();j++){ // loop over hits in sub-event
					TTrificHit* hitb=TimeCorrelatedSingles[k][j];
					if(abs(hit->GetTime()-hitb->GetTime())>100)break; // If time to different move to next sub-event
					if(j+1==TimeCorrelatedSingles[k].size()){// If passed time check with all hits in sub-events, add hit to sub-evet
						UnUsed=false;
						TimeCorrelatedSingles[k].push_back(hit);
                        
                        if(tE>dEMax[k])dEMax[k]=tE;
                        if(d>Range[k])Range[k]=d;
                        Contin[k]*=d;
                        if(d==2)dE2[k]=tE;
                        if(d==4)dE4[k]=tE;
						break;
					}
				}
				if(!UnUsed) break;
			}
			
			// If hit not added to existing sub-event, create new sub event
			if(UnUsed){
				TimeCorrelatedSingles.push_back(vector<TTrificHit*>{hit});
                
                dEMax.push_back(tE);
                Range.push_back(d);
                Contin.push_back(d);
                if(d==2)dE2.push_back(tE);else dE2.push_back(0);
                if(d==4)dE4.push_back(tE);else dE4.push_back(0);
			}
		}
		
		// Multiplicity checks
        for(int i=0;i<24;i++){
			SingleMult->Fill(i,Ncount[i]);
			if(Ncount[i]){
				for(int j=0;j<24;j++){
					if(Ncount[j]){
						TrigCount->Fill(i,j);
					}
				}
			}
		}
		
		
		
        if(T6){ 
            for(int h=0;h<trific->GetSingMultiplicity();h++){
                TTrificHit* hit=trific->GetTrificSingHit(h);
                int d=hit->GetDetector();
                if(d==6)continue;
                dT6->Fill(T6-hit->GetTime(),d);
            }
        }
        
		
		// Loop over X grids
        for(int h=0;h<trific->GetXMultiplicity();h++){
            TTrificHit* hit=trific->GetTrificXHit(h);
            double tE=hit->GetEnergy();
            int s=hit->GetSegment();
            XSingles->Fill(s,tE);
			XkValue->Fill(s,log10(hit->GetKValue()));
            EX->Fill(s,tE);
            
            if(T6){
                dTEX6->Fill(T6-hit->GetTime(),tE);
            }
		}
		
		
		// Loop over Y grids
        for(int h=0;h<trific->GetYMultiplicity();h++){
            TTrificHit* hit=trific->GetTrificYHit(h);
            double tE=hit->GetEnergy();
            int s=hit->GetSegment();
            YSingles->Fill(s,tE);
			YkValue->Fill(s,log10(hit->GetKValue()));
            
            EY->Fill(s,tE);
            YN->Fill(s);
            if(T6){
                YN6->Fill(s);
                dTimeY6->Fill(T6-hit->GetTime());
                dTEY6->Fill(T6-hit->GetTime(),tE);
            }
            
			for(int k=0;k<trific->GetXMultiplicity();k++){
				TTrificHit* hitb=trific->GetTrificXHit(k);
				XYMAPRAW->Fill(hitb->GetSegment(),hit->GetSegment());
				
				if(trific->GetYMultiplicity()==1&&trific->GetXMultiplicity()==1&&Ncount[6]>0){
					XYmap_single->Fill(hitb->GetSegment(),hit->GetSegment());
				}
				
			}
		}
		
		XYMult->Fill(trific->GetXMultiplicity(),trific->GetYMultiplicity());
		
        for(int h=0;h<trific->GetSingMultiplicity();h++){
            TTrificHit* hit=trific->GetTrificSingHit(h);
            int d=hit->GetDetector();
            double tE=hit->GetEnergy();
            if(tE10>100){
                Singles10->Fill(d,tE);
            }
            
            if(tE12>100){
                Singles12->Fill(d,tE);
            }
        }

		

		for(unsigned int k=0;k<TimeCorrelatedSingles.size();k++){ //loop over existing sub-events

            
            // If continuity 
            double R=Range[k];
            if(R>23)continue;
            if((Rfact[R]==Contin[k])||(Rfact[R]==Contin[k]*2)||(Rfact[R]==Contin[k]*9)||(Rfact[R]==Contin[k]*18)){
                 
                ID_RangeZ->Fill(R,dEMax[k]);
                ID_Range4->Fill(R,dE4[k]);
                ID_Range2->Fill(R,dE2[k]);
                ID_Z2->Fill(dEMax[k],dE2[k]);
                ID_Z4->Fill(dEMax[k],dE4[k]);
                
                for(unsigned int h=0;h<TimeCorrelatedSingles[k].size();h++){ // loop over hits in sub-event
                    TTrificHit* hit=TimeCorrelatedSingles[k][h];
                    int d=hit->GetDetector();
                    double tE=hit->GetEnergy();
                    ContinuitySingles->Fill(d,tE);
                }

            }
           
            
            
			for(unsigned int h=0;h<TimeCorrelatedSingles[k].size();h++){ // loop over hits in sub-event
				TTrificHit* hit=TimeCorrelatedSingles[k][h];
				int d=hit->GetDetector();
				double tE=hit->GetEnergy();
				
				
				for(unsigned int j=h+1;j<TimeCorrelatedSingles[k].size();j++){ // loop over hits in sub-event
					TTrificHit* hitb=TimeCorrelatedSingles[k][j];
					dTGated->Fill(hitb->GetTime()-hit->GetTime());
				}
			}
		}
    }
    
	out.Write();
	out.Close();
}

	

