{
    //TChannel::ReadCalFile("trificCal.cal");
    TChannel::ReadCalFile("CalibrationFile.cal");

    
    // Set false after parameter are tuned to data
    bool DoTestDraw=true;

    std::cout << "\nPlease enter the run number: ";
    string input;
    std::cin >> input;

    std::string stringFile = "/tig/pterodon_data3/S1893/AnalysisTrees/analysis"+input+"_000.root";
    //std::string stringFile = "/tig/pterodon_data3/S1893/AnalysisTrees/analysis"+input+"_000.root";
    const char* inFile = stringFile.c_str();

    
    // Loading analysis tree etc (this is all done for you in grsiproof)
 ///////// FILL IN A FEW MANUAL ANALYSIS TREE FILE NAMES ///////// 
    TChain* DataChain = new TChain("AnalysisTree","AnalysisTree");
//     DataChain->Add("analysis51989_000.root");
    //DataChain->Add("/tig/pterodon_data3/S1893/testAnalysisTrees/analysis51994_000.root");
	


	DataChain->Add(inFile);    



    // Set address of branch
    TTrific *trific=0;
    if(DataChain->FindBranch("TTrific"))DataChain->SetBranchAddress("TTrific",&trific);else return;

    // Load things and set things for class statics
    TTrific::SetSharc();
//     TVector3 gridnormal=TTrific::fNormalGridVec;
    double fAngle = (60./180.)*TMath::Pi();
    TVector3 gridnormal=TVector3(0,-TMath::Cos(fAngle),TMath::Sin(fAngle));
    
    double fInitialSpacingCart=TTrific::fInitialSpacingCart;
    double fSpacingCart=TTrific::fSpacingCart;

 ///////// SHOULD THIS BE ZERO? I COULDNT REMEMBER ///////// 
    int FirstGridN=1;

    
    // Create the TF1 fit functions we're going to use in the loop
 ///////// THESE PARAMETER ESTIMATES NEED TUNING TO THE DATA SET ///////// 
    TF1 tail("tail","pol1",0,20);
    tail.SetParameters(14000,1000);
    tail.SetParNames("offset","gradient");
    
    TF1 peak("peak","[1]-(x<[0])*abs([2])*pow(x-[0],2)-(x>=[0])*abs([3])*pow(x-[0],2)");
    peak.SetParameters(32,20000,1000,10000);
    peak.SetParNames("peakX","peakY","lowsig","highsig");

    TF1 newbragg("newbragg","[0]*exp(([1]*[1]*[2]*[2]*0.5)+([1]*(x-[3])))*TMath::Erfc(([1]*[2]/sqrt(2))-(([3]-x)/([2]*sqrt(2))))",0,50);
    newbragg.SetLineColor(3);

    // Quick example histogram
    TH2F* rangeZ=new TH2F("rangeZ","Range vs dE/dx peak;Range [cm];de/dx peak [arb.]",200,20,40,5000,0,10000);
    TH2F* rangeZextra=new TH2F("rangeZextra","Range vs dE/dx peak Full Fit;Range [cm];de/dx peak [arb.]",200,20,40,5000,0,10000);
    
    // Loop over events
    Int_t nentries = DataChain->GetEntries();
    for(int jentry = 0; jentry<nentries*0.1; jentry++){
        DataChain->GetEntry(jentry);

        // Get Position, skip event if no good
        TVector3 particle=trific->GetPosition();
        if(particle.Z()<=0)continue;
        
        TGraph tmpbragg;
        double maxE=0;
        double maxX=0;
        double ZratioX=1./gridnormal.Dot(particle.Unit());

        // Loop over hits
        for(int h=0;h<trific->GetSingMultiplicity();h++){
            TTrificHit* hit=trific->GetTrificSingHit(h);
            double tE=hit->GetEnergy();
            
            // Skip "bad" hits
            if(tE<0.1)continue;
            if(hit->GetDetector()==FirstGridN)continue;
            
            // Calculated de/dx and "corrected Z"
            double basicZ=fInitialSpacingCart+fSpacingCart*(hit->GetDetector()-FirstGridN);
            basicZ*=0.1;//My function estimates are in cm not mm
            double dE=tE/ZratioX;
            double x=basicZ*ZratioX;
            tmpbragg.SetPoint(tmpbragg.GetN(),x,dE);
            
            if(dE>maxE){
                maxE=dE;
                maxX=x;
            }
        }
 ///////// I HAVE NOT IMPLEMENTED IT HERE BUT WE SHOULD ALSO CREATE A DATA POINT FOR THE X AND Y GRIDS BY SUMMING THEIR ENERGIES ///////// 
        
        if(tmpbragg.GetN()<4)continue;
        
        // Fit a very basic function to the bragg peak to get (de/dx)_peak and range
 ///////// THESE PARAMETER AND RANGE ESTIMATES NEED TUNING TO THE DATA SET ///////// 
        peak.SetParameters(maxX,maxE,1000,10000);
        peak.SetRange(maxX-5,40);
        tmpbragg.Fit(&peak,"QNR");   
    
        // Extract estimated values from fit
        double dedxpeak=peak.GetParameter(1);
        double range=peak.GetParameter(0)+sqrt(abs(dedxpeak/peak.GetParameter(3))); 
        
        // Fit a very basic function to the bragg peak to get (de/dx)_peak and range
 ///////// THESE PARAMETER AND RANGE ESTIMATES NEED TUNING TO THE DATA SET ///////// 
//         tail.SetParameters(14000,1000);
//         tmpbragg.Fit(&tail,"QNR");      
//         double dE0=tail.Eval(5);
        
        rangeZ->Fill(range,dedxpeak);

        // Fit the more complicated full curve fit
        double xmax=peak.GetParameter(0);
        double sig=range-xmax;
        if(sig<1)sig=1;
        newbragg.SetParameter(0,dedxpeak*0.5);
        newbragg.SetParLimits(0,dedxpeak*0.4,dedxpeak*0.6);
        newbragg.SetParameter(1,0.02);
        newbragg.SetParameter(2,sig);
        newbragg.SetParLimits(2,sig*0.1,sig*2);
        newbragg.SetParameter(3,xmax+sig*0.5);
        newbragg.SetParLimits(3,xmax-sig,xmax+sig*2);
        tmpbragg.Fit(&newbragg,"QNR");
         
        double dE0=newbragg.Eval(5);
        range=newbragg.GetParameter(3)+newbragg.GetParameter(2)*2;// Approx
        dedxpeak=newbragg.Eval(newbragg.GetParameter(3)-newbragg.GetParameter(2)*2.3);// Approx
        
        rangeZextra->Fill(range,dedxpeak);
        
        
        // Draw the fit results to check things are working
        if(DoTestDraw){
            TCanvas C1;
            gPad->Update();
            tmpbragg.Fit(&tail,"+R");
            tmpbragg.Fit(&peak,"+R");
            tmpbragg.Fit(&newbragg,"+R");
            tmpbragg.SetMarkerStyle(20);
            tmpbragg.Draw("ap");        
            
            cout<<endl<<"dE0 = "<<dE0;
            cout<<endl<<"Z = "<<dedxpeak;
            cout<<endl<<"range = "<<range<<endl;
            
            C1.Modified();
            C1.Update();
            C1.WaitPrimitive();
        }
        
        if(jentry%1000 == 0){
            cout << setiosflags(ios::fixed) << std::setprecision(2) << 100.0*jentry/nentries << " % complete.\r" << flush;
        }
    }
    
    // Draw the ID plots
    TCanvas* C2=new TCanvas();
    C2->Divide(2);
    C2->cd(1);
    rangeZ->Draw("col");
    C2->cd(2);
    rangeZextra->Draw("col");
}

	

