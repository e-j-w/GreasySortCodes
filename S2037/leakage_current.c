//Leakage Current Monitoring
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "time.h"

void leakage_current() {
  float S3uA, SSB1uA, SSB2uA, FPuA;
  time_t now = time(0);
  //cout << "S3 Leakage current \u03BCA = ";
  //cin >> S3uA;
  cout << "SSB Left Leakage current \u03BCA = ";
  cin >> SSB1uA;
  cout << "SSB Right Leakage current \u03BCA = ";
  cin >> SSB2uA;
  cout << "Focal Plane Leakage current \u03BCA = ";
  cin >> FPuA;

  //string resultfile = "Monitoring/S3LeakageCurrent.dat";
  //ifstream test(resultfile.c_str());
  //ofstream s3file;
  //if(test.good()) {
    //s3file.open(resultfile.c_str(),ios::app);
  //}
  //else {
    //s3file.open(resultfile.c_str());
  //}
  //s3file << now << "\t" << S3uA << endl;
  //s3file.close();

  string resultfile = "Monitoring/SSBLeftLeakageCurrent.dat";
  ifstream test2(resultfile.c_str());
  ofstream ssb1file;
  if(test2.good()) {
    ssb1file.open(resultfile.c_str(),ios::app);
  }
  else {
    ssb1file.open(resultfile.c_str());
  }
  ssb1file << now << "\t" << SSB1uA << endl;
  ssb1file.close();

  resultfile = "Monitoring/SSBRightLeakageCurrent.dat";
  ifstream test3(resultfile.c_str());
  ofstream ssb2file;
  if(test3.good()) {
    ssb2file.open(resultfile.c_str(),ios::app);
  }
  else {
    ssb2file.open(resultfile.c_str());
  }
  ssb2file << now << "\t" << SSB2uA << endl;
  ssb2file.close();

  resultfile = "Monitoring/FocalPlaneLeakageCurrent.dat";
  ifstream test4(resultfile.c_str());
  ofstream fpfile;
  if(test4.good()) {
    fpfile.open(resultfile.c_str(),ios::app);
  }
  else {
    fpfile.open(resultfile.c_str());
  }
  fpfile << now << "\t" << FPuA << endl;
  fpfile.close();

  //TCanvas *c1 = new TCanvas();
  //c1->Draw();
  //c1->cd();
  //TGraph *gr = new TGraph("Monitoring/S3LeakageCurrent.dat");
  //gr->Draw("A*");
  //gr->GetXaxis()->SetTimeDisplay(1);
  //gr->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  //gr->GetYaxis()->SetTitle("Leakage Current (#muA)");

  TCanvas *c2 = new TCanvas();
  c2->Draw();
  c2->cd();
  TGraph *gr2 = new TGraph("Monitoring/SSBLeftLeakageCurrent.dat");
  gr2->Draw("A*");
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr2->GetYaxis()->SetTitle("Leakage Current (#muA)");

  TCanvas *c3 = new TCanvas();
  c3->Draw();
  c3->cd();
  TGraph *gr3 = new TGraph("Monitoring/SSBRightLeakageCurrent.dat");
  gr3->Draw("A*");
  gr3->GetXaxis()->SetTimeDisplay(1);
  gr3->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr3->GetYaxis()->SetTitle("Leakage Current (#muA)");

  TCanvas *c4 = new TCanvas();
  c4->Draw();
  c4->cd();
  TGraph *gr4 = new TGraph("Monitoring/FocalPlaneLeakageCurrent.dat");
  gr4->Draw("A*");
  gr4->GetXaxis()->SetTimeDisplay(1);
  gr4->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr4->GetYaxis()->SetTitle("Leakage Current (#muA)");

}
