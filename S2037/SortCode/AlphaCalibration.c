//g++ AlphaCalibration.c -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o Alphacal
#include <iostream> 
#include <iomanip> 
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdio.h> 
#include <stdlib.h> 
#include <TCanvas.h>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TTip.h"
#include "TTigress.h"
#include "TRF.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

void HistDraw(char const* infile, char const* calfile, char const* outfile) {
    
    Double_t si_nBins = 4096;
    Double_t si_min = 0;
    Double_t si_max = 1024*16;

    TList * list = new TList;
    //TH2F *energy_channelN = new TH2F("energy_channelN","Energy vs Channel Number",4096,0,4096,20,0,20);list->Add(energy_channelN);
    TH1F *myhist[5000];
    char histname[20];
    for (int iii=0;iii<2000; iii++) {
      sprintf(histname, "Energy_ch%1.1i", iii); //Channel%d
      myhist[iii] = new TH1F(histname, Form("myhist d%1.1i", iii), si_nBins, si_min, si_max);
    }
    

    TFile * inputfile = new TFile(infile, "READ");
    if (!inputfile->IsOpen()) {
      printf("Opening file failed, aborting\n");
      return;
    }
    TChain * FragmentTree = (TChain * ) inputfile->Get("FragmentTree");
    printf("%i tree files, details:\n", FragmentTree->GetNtrees());
    FragmentTree->ls();
    TTree * tree = (TTree * ) FragmentTree->GetTree();
    Int_t nentries = FragmentTree->GetEntries();
    TFragment * frag = 0;
    FragmentTree->SetBranchAddress("TFragment", & frag);
    printf("Reading calibration file: %s\n", calfile);
    TChannel::ReadCalFile(calfile);
    int npeaks = 20;

    printf("Begin sort\n");
    int one;
    for (int jentry = 0; jentry < (nentries); jentry++) {
      tree->GetEntry(jentry);
      if (frag->GetChannelNumber() >= 0) myhist[frag->GetChannelNumber()]->Fill(frag->GetEnergy());
      if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Si Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
    }

    for (int iii = 0; iii < 2000; iii++) {
      //for(int jjj=0; jjj<200; jjj++)myhist[iii]->SetBinContent(jjj, 0);
      if (myhist[iii]->Integral(0, si_nBins) > 10) list->Add(myhist[iii]);
    }

        cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

        cout << "Histograms written to Hist.root, sorting complete" << endl;

        TFile * myfile = new TFile(outfile, "RECREATE");
        myfile->cd();
        list->Write();
	myfile->Write();
        myfile->Close();
}

Double_t tripleAlphaSpectrum(Double_t *xx, Double_t *p){
// Simulate a triple alpha spectrum for comparison to a histogram
// Energy calibration and width of peaks are parameters given up to quadratic linear
// Parameters:
// p0 -- Normalization factor for 239Pu group
// p1 -- Normalization factor for 241Am group
// p2 -- Normalization factor for 244Cm group
// p3 -- FWHM of peaks in keV
// p4 -- Constant offset in counts*
// p5 -- Width of Pu peaks relative to Cm
// p6-- Width of Am peaks relative to Cm
// p7 -- Offset
// p8-- Linear gain
 Double_t E = p[7]+p[8]*xx[0];
 Double_t sigmaCu = p[3]/2.35;
 Double_t sigmaPu = sigmaCu*p[5];
 Double_t sigmaAm = sigmaCu*p[6];
// Now calculate spectrum height based on energies and proportionalities
 Double_t f = p[0] * ( 0.7077 * TMath::Gaus(E,5156.59,sigmaPu)
          +0.1711 * TMath::Gaus(E,5144.30,sigmaPu)
          +0.1194 * TMath::Gaus(E,5105.80,sigmaPu) )
        +p[1] * ( 0.0036 * TMath::Gaus(E,5544.5,sigmaAm)
           +0.0166 * TMath::Gaus(E,5388,sigmaAm)
           +0.848 * TMath::Gaus(E,5485.56,sigmaAm)
           +0.131 * TMath::Gaus(E,5442.8,sigmaAm) )
        +p[2] * ( 0.769 * TMath::Gaus(E,5804.77,sigmaCu)
           +0.231 * TMath::Gaus(E,5762.16,sigmaCu ) )
        +p[4];
  return f;
}
TF1* tasf(TH1* h, const char* name = "tas"){
// Create a TF1 using the function above and set some intial parameters
 TF1* F = new TF1(name,tripleAlphaSpectrum,0,10000,9); //this sets the limits of the fit in channels (see also line 93)
 F->SetParNames("Pu","Am","Cm","fwhmCm","a","Pu_n","Am_n","offset","gain");
 F->SetNpx(7000);
 F->SetParLimits(0,10,500000);
 F->SetParLimits(1,10,500000);
 F->SetParLimits(2,10,500000);
 F->SetParLimits(3,20,200);
 F->SetParLimits(5,0.8,1.3);
 F->SetParLimits(6,0.8,1.2);
// Search for the peaks in the charge spectrum to initialize gain and offset
 vector< Double_t > peaks;
 int npeaks = 20;
 TSpectrum * s = new TSpectrum(2*npeaks);
 Int_t nfound = s->Search(h, 2, "", 0.13);
 if (nfound>=3){
    Double_t * xpeaks = s->GetPositionX();  
    for (int oo=0; oo<=nfound; oo++){
       if (xpeaks[oo]>400 and xpeaks[oo]<10000){
          peaks.push_back(xpeaks[oo]); 
          cout << "peak: " << xpeaks[oo] << endl;  
      }
    }
    auto min = *std::min_element( peaks.begin(), peaks.end());
    auto max = *std::max_element( peaks.begin(), peaks.end());
    Double_t m = (5804.77-5156.59)/(max-min);
    Double_t b = 5804.77-m*max; 
//Initializing Parameters
    F->SetParameters(50,50,50,25,0,1,1,b,m);
    F->SetParLimits(7,b-100,b+100);
    F->SetParLimits(8,m-0.4,m+0.4);
 }
 return F;
}
TH1* getHistFromFile(const char* fn, const char* objpath){
 TFile* F = TFile::Open(fn);
 TFolder* Fo;
 TH1* H;
 F->GetObject(objpath,H);
 H->SetTitle(fn);
 H->Draw("hist && E");
 return H;
}
TF1* zoomThree(TH1* h, TF1* F){
// Function to zoom h around the peaks for fitting
 Double_t g = F->GetParameter("gain");
 Double_t o = F->GetParameter("offset");
 
 h->GetXaxis()->SetRange((800-o)/g,(2500-o)/g); //this sets the region to zoom the plot into
 h->Draw();
 return F;
}
void Peaks(Double_t ratio, Double_t constant, Double_t mean, Double_t sigma){
// Function to draw the individual peaks of the main fit
   TF1 *myfit = new TF1("myfit","gaus",0,10000); //this sets the limits of the fit in channels (see also line 41)
   myfit->SetParameters(ratio * constant , mean , sigma);
   myfit->SetLineColor(3);
   myfit->Draw("SAME");
}

// Find resolutions of Cm for the working channels using the functions above
// Only channels start to end 
void m( int start, int end){
   gROOT->SetBatch(kTRUE);
   vector< double > workingchannels;
   vector< double > warr;
   vector< double > werrarr;
   vector< double > rchiarr;
   vector< double > OFFSET;
   vector< double > GAIN;
   TList * list = new TList;
   TCanvas * C[end-start];
   //string histfile;
   //cout << "\nEnter the name of the histogram file.\n";
   //getline(cin, histfile);
   for (int iii=start; iii<=end; iii++){
// Create canvases that will be output to TAS.root
      C[iii-start] = new TCanvas(Form("c%1.1i",iii),Form("c%1.1i",iii),30,113,
      800,600);
      C[iii-start]->Range(0,-10,60,600);
// Gather, draw, zoom around, and fit data. To display fit parameters, remove Q
// Open the histograms by objpath (ie. Energy_ch##) within fmc32.root (fn)
// For analyzing SHARC data, (ie. Channel##) with oldsharc.root or newsharc.root
      //string hist;
     // cout << "\nEnter the histogram file name.\n";
      //getline(cin, hist);
      TH1* h = getHistFromFile("Hist.root",Form("Energy_ch%1.1i",iii));  
      cout << iii << " :";   
      TF1* f=tasf(h); 
      f = zoomThree(h,f); f->Draw("SAME"); f = zoomThree(h,f); 
      h->Fit("tas","LQ"); h->Fit("tas","LQ");
// Get parameters, FWHM, and draw the individual peaks    
      Double_t g = f->GetParameter("gain");
      Double_t o = f->GetParameter("offset");     
      Double_t fwhmCm = f->GetParameter("fwhmCm");
      Double_t fwhmCuerr = f->GetParError(3);
      Double_t rchi2 = f->GetChisquare()/f->GetNDF();
      Double_t sigmaCm = fwhmCm/2.35;
      Double_t Pu_n = f->GetParameter("Pu_n"); Double_t sigmaPu = (Pu_n*sigmaCm);
      Double_t Am_n = f->GetParameter("Am_n"); Double_t sigmaAm = (Am_n*sigmaCm);
      Double_t Cm = f->GetParameter("Cm");
      Double_t Pu = f->GetParameter("Pu");
      Double_t Am = f->GetParameter("Am");
      Double_t ratios[9]={0.7077,0.1711,0.1194,0.0036,0.0166,0.848,0.131,0.769,
      0.231};
      Double_t constants[9]={Pu,Pu,Pu,Am,Am,Am,Am,Cm,Cm};           
      Double_t means[9]={(5156.59-o)/g,(5144.30-o)/g,(5105.80-o)/g,(5544.5-o)/g,
      (5388-o)/g,(5485.56-o)/g,(5442.8-o)/g,(5804.77-o)/g,(5762.16-o)/g};
      Double_t sigmas[9]={sigmaPu/g,sigmaPu/g,sigmaPu/g,sigmaAm/g,sigmaAm/g,
      sigmaAm/g,sigmaAm/g,sigmaCm/g,sigmaCm/g};
      for (int ii=7; ii<8; ii++){
         Peaks(ratios[ii],constants[ii],means[ii],sigmas[ii]);
      }
// Add canvas to list      
      list->Add(C[iii-start]);  
// Add the FWHM and error of Curium peak
      workingchannels.push_back(iii);
      warr.push_back(fwhmCm);
      werrarr.push_back(fwhmCuerr);
      rchiarr.push_back(rchi2);
      OFFSET.push_back(o);
      GAIN.push_back(g);
   }
// Output the widths of Curium for each channel and its absolute error
   cout << "CHANNEL" << "\t" << "FWHM (keV)" << "\t\t" << "ERR" << "\t\t" 
   << "REDUCED CHISQR" << endl;
// Output gains and offsets to a text file
   ofstream file;
   file.open ("Calibration.txt");
   file << "float GAIN[" << end - start + 1 << "] = {" << GAIN[0];
   for (int j=0; j<warr.size(); j++){    
      cout << workingchannels[j] << ",\t" << warr[j] << ",\t\t" << werrarr[j] 
      << ",\t" << rchiarr[j] << endl;
      if (j > 0) {
         file << ", " << GAIN[j];
      }
   } 
   file << "};\n" << "float OFFSET[" << end - start + 1 << "] = {" << OFFSET[0];
   for (int j=0; j<warr.size(); j++){    
      if (j > 0) {
         file << ", " << OFFSET[j];
      }
   } 
   file << "};"; 
   file.close();
   cout << "Writing gains and offsets to " << "Calibration.txt" << endl;
   
// Write canvases to TAS.root
   cout << "Writing histograms to " << "TAS.root" << endl;
   TFile * myfile = new TFile("TAS.root", "RECREATE");
   myfile->cd();
   list->Write();
   myfile->Close();
}

int main(int argc, char **argv){
	
	char const *ffile;
	char const *calfile;
	char const *outfile;
	char const *ChMin;
        char const *ChMax;
	int ChMin_int;
	int ChMax_int;

	printf("Starting sortcode\n");

  	std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  	if(grsi_path.length() > 0) {
	  grsi_path += "/";
 	}
  	// Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
	grsi_path += ".grsirc";
  	gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

	if(argc == 1)
	{
		cout << "Insufficient arguments, provide argument tree files" << endl;
		return 0;
	}

	else if(argc == 2)
	{
		ffile = argv[1];
		calfile = "/CalibrationFile.cal";
		ChMin = "0";
	        ChMax = "55";
	}
	else if(argc == 3)
	{
		ffile   = argv[1];
		calfile = argv[2];
		ChMin = "0";
		ChMax = "55";
	}
	else if(argc == 4)
	{
		ffile   = argv[1];
		calfile = argv[2];
		ChMin = argv[3];
		ChMax = "55";
	}
	else if(argc == 5)
	{
		ffile   = argv[1];
		calfile = argv[2];
		ChMin = argv[3];
		ChMax = argv[4];
	}
	else if(argc > 5)
	{
		printf("Too many arguments\n");
		return 0;
	}

	outfile = "Hist.root";

	printf("Input file:%s\nCalibration file: %s\nOutput File: %s\nStarting Channel: %s\nEnding Channel: %s\n",ffile,calfile,outfile,ChMin,ChMax);
  	TParserLibrary::Get()->Load();
	
	//Convert the channel minimum and maximum characters to integers so they can be used as input in m()
	ChMin_int = atoi(ChMin);
	ChMax_int = atoi(ChMax);

	HistDraw(ffile,calfile, outfile);
	m(ChMin_int,ChMax_int);
	return 0;
}

