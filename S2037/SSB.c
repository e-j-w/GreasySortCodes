//SSB Monitoring Script
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"

void SSB() {

  int z = 52880; // First Run on target
  TH2F *h1[100];
  TH2F *h2[100];
  TH2F *h3[100];
  double rnum;
  double rnumer;
  double SSB1,SSB2,REC1;
  double SSB1er,SSB2er,REC1er;
  double temp_current;
  int j = 0;
  int k = 0;
  int bad_runs[] = {52816}; //JW: there needs to be at least one run here to avoid runtime errors

  Int_t rNum[300]; 
  Double_t bC[300];
  string rfile = "Monitoring/beam_currents.dat";
  ifstream fp;
  fp.open(rfile.c_str());
  int iii =0;
  if(!fp.is_open())cout << "File Not Found" << endl;
  while(fp.good()) {
    fp >> rNum[iii] >> bC[iii];
    //printf("run %i current %f\n",rNum[iii],bC[iii]);
    iii++; 
  }
  fp.close();

  ifstream test(rfile.c_str());
  ofstream fpout;
  if(test.good()) {
    fpout.open(rfile.c_str(),ios::app);
  }
  else {
    fpout.open(rfile.c_str());
  }

  string ssb1 = "Monitoring/SSBLeft.dat";
  string ssb2 = "Monitoring/SSBRight.dat";
  string rec1 = "Monitoring/Recoils.dat";
  ifstream test2(ssb1.c_str());
  ofstream ssb1file;
  if(test2.good()) {
    ssb1file.open(ssb1.c_str(),ios::app);
  }
  else {
    ssb1file.open(ssb1.c_str());
  }
  ifstream test3(ssb2.c_str());
  ofstream ssb2file;
  if(test3.good()) {
    ssb2file.open(ssb2.c_str(),ios::app);
  }
  else {
    ssb2file.open(ssb2.c_str());
  }
  ifstream test4(rec1.c_str());
  ofstream rec1file;
  if(test4.good()) {
    rec1file.open(rec1.c_str(),ios::app);
  }
  else {
    rec1file.open(rec1.c_str());
  }

  for(int i=0;i<300;i++) {
    int * ap = find(begin(bad_runs), end(bad_runs), z+i); 
    if (ap != end(bad_runs)) continue;
    if(gSystem->AccessPathName(Form("HistFiles/Hist_%i_000.root",(z+i)))){
      continue;
    } else { 
      int *q = find(begin(rNum), end(rNum), z+i); 
//cout << "rNum" << rNum << " q " << q << " end(rNum) " << end(rNum) << endl;
      if (q == end(rNum)) {
        cout << "Beam Current (run " << (z+i) << ") (nA) = ";
        cin >> temp_current; 
        TFile *infile;
        infile = new TFile(Form("HistFiles/Hist_%i_000.root",(z+i)));
        infile->cd();
        h1[j] = (TH2F*)infile->Get("EMMA/SSB_1_Versus_Time");
        h2[j] = (TH2F*)infile->Get("EMMA/SSB_2_Versus_Time");
        h3[j] = (TH2F*)infile->Get("EMMA/icSumVSi");
        rnum = z+i;
        rnumer = 0;
        SSB1 = h1[j]->Integral(100,400,0,15000)/temp_current; //y axis values need to be divided by 6 to convert from keV to channels, due to changes made by JW to binning of these histograms in SortCode.h, x axis still has 1:1 mapping
        SSB1er = pow(h1[j]->Integral(100,400,0,160000),0.5)/h1[j]->Integral(100,400,0,160000) * SSB1;
        SSB2 = h2[j]->Integral(100,400,0,15000)/temp_current; //y axis values need to be divided by 6 to convert from keV to channels, due to changes made by JW to binning of these histograms in SortCode.h, x axis still has 1:1 mapping
        SSB2er = pow(h2[j]->Integral(100,400,0,160000),0.5)/h2[j]->Integral(100,400,0,160000) * SSB2;
        //REC1 = h3[j]->Integral(1300,2200,1800,2300)/temp_current; //y axis values need to be divided by 6 to convert from keV to channels, due to changes made by JW to binning of these histograms in SortCode.h, x axis still has 1:1 mapping
       // REC1er = pow(h3[j]->Integral(1200,2200,1800,2200),0.5)/h3[j]->Integral(1200,2200,1800,2200) * SSB1;
        j++;
        infile->Close();
	fpout << rnum << "\t" << temp_current << endl;
	ssb1file << rnum << "\t" << SSB1 << "\t" << SSB1er << endl;
        ssb2file << rnum << "\t" << SSB2 << "\t" << SSB2er << endl;
        //rec1file << rnum << "\t" << REC1 << "\t" << REC1er << endl;
      }
      else continue;
    }
  }
  fpout.close();
  ssb1file.close();
  ssb2file.close();
  rec1file.close();
  TGraphErrors *ssb1g = new TGraphErrors("Monitoring/SSBLeft.dat","%lg %lg %lg");
  ssb1g->SetLineColor(kBlack);
  ssb1g->SetMarkerColor(kBlack);
  ssb1g->SetTitle("Left SSB");
  ssb1g->SetFillColor(0);
  TGraphErrors *ssb2g = new TGraphErrors("Monitoring/SSBRight.dat","%lg %lg %lg");
  ssb2g->SetLineColor(kRed);
  ssb2g->SetMarkerColor(kRed);
  ssb2g->SetTitle("Right SSB");
  ssb2g->SetFillColor(0);
/*  TGraphErrors *rec1g = new TGraphErrors("Monitoring/Nov2020/Recoils.dat","%lg %lg %lg");
  rec1g->SetLineColor(kBlue);
  rec1g->SetMarkerColor(kBlue);
  rec1g->SetTitle("Recoils");
  rec1g->SetFillColor(0);
*/  TCanvas *c1 = new TCanvas("c1","",600,400);
  c1->cd();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);

  TMultiGraph *gr = new TMultiGraph();
  gr->Add(ssb1g);
  gr->Add(ssb2g);
  gr->Draw("A*");

// Adding vertical lines to show when target has been changed
  int target_changes[] = {z};

  TLine* line = new TLine();
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  for(int tc_iter=0; tc_iter<sizeof(target_changes)/sizeof(target_changes[0]); tc_iter++){
    double_t x_pos = target_changes[tc_iter]+0.5;
    line->DrawLine(x_pos, gr->GetYaxis()->GetXmin(), x_pos, 12.);
  }

  TLegend *leg = new TLegend(0.65,0.65,0.85,0.88);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(ssb1g);
  leg->AddEntry(ssb2g);
  leg->Draw();
/*
  //make a 2nd plot
  c1 = new TCanvas("c1","",600,400);
  c1->cd();
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);

  gr = new TMultiGraph();
  gr->Add(rec1g);
  gr->Draw("A*");

  leg = new TLegend(0.65,0.65,0.85,0.88);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.045);
  leg->AddEntry(rec1g);
  leg->Draw();

  line = new TLine();
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  for(int tc_iter=0; tc_iter<sizeof(target_changes)/sizeof(target_changes[0]); tc_iter++){
    double_t x_pos = target_changes[tc_iter]+0.5;
    line->DrawLine(x_pos, gr->GetYaxis()->GetXmin(), x_pos, 12.);
  }
*/
}
