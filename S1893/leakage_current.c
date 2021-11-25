//Leakage Current Monitoring
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "time.h"

void leakage_current() {
  float S3uA, SSB1uA, SSB2uA, FPuA;
  float shb05, shb06, shb07, shb08, shb09, shb10, shb11, shb12, shq13, shq14, shq15, shq16, pad1, pad2, pad3, pad4, trifPos, trifNeg; //all in uA
  time_t now = time(0);
  cout << "SHB05 Leakage current \u03BCA = ";
  cin >> shb05;
  cout << "SHB06 Leakage current \u03BCA = ";
  cin >> shb06;
  cout << "SHB07 Leakage current \u03BCA = ";
  cin >> shb07;
  cout << "SHB08 Leakage current \u03BCA = ";
  cin >> shb08;
  cout << "SHB09 Leakage current \u03BCA = ";
  cin >> shb09;
  cout << "SHB10 Leakage current \u03BCA = ";
  cin >> shb10;
  cout << "SHB11 Leakage current \u03BCA = ";
  cin >> shb11;
  cout << "SHB12 Leakage current \u03BCA = ";
  cin >> shb12;

  cout << "SHQ13 Leakage current \u03BCA = ";
  cin >> shq13;
  cout << "SHQ14 Leakage current \u03BCA = ";
  cin >> shq14;
  cout << "SHQ15 Leakage current \u03BCA = ";
  cin >> shq15;
  cout << "SHQ16 Leakage current \u03BCA = ";
  cin >> shq16;

  cout << "Pads 1 Leakage current \u03BCA = ";
  cin >> pad1;
  cout << "Pads 2 Leakage current \u03BCA = ";
  cin >> pad2;
  cout << "Pads 3 Leakage current \u03BCA = ";
  cin >> pad3;
  cout << "Pads 4 Leakage current \u03BCA = ";
  cin >> pad4;

  cout << "TRIFIC Positive Bias Leakage current \u03BCA = ";
  cin >> trifPos;
  cout << "TRIFIC Negative Bias Leakage current \u03BCA = ";
  cin >> trifNeg;

  

  string resultfile = "Monitoring/SHB05LeakageCurrent.dat";
  ifstream test(resultfile.c_str());
  ofstream shb05file;
  if(test.good()) {
    shb05file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb05file.open(resultfile.c_str());
  }
  shb05file << now << "\t" << shb05 << endl;
  shb05file.close();

  resultfile = "Monitoring/SHB06LeakageCurrent.dat";
  ifstream test2(resultfile.c_str());
  ofstream shb06file;
  if(test2.good()) {
    shb06file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb06file.open(resultfile.c_str());
  }
  shb06file << now << "\t" << shb06 << endl;
  shb06file.close();

  resultfile = "Monitoring/SHB07LeakageCurrent.dat";
  ifstream test3(resultfile.c_str());
  ofstream shb07file;
  if(test3.good()) {
    shb07file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb07file.open(resultfile.c_str());
  }
  shb07file << now << "\t" << shb07 << endl;
  shb07file.close();

  
  resultfile = "Monitoring/SHB08LeakageCurrent.dat";
  ifstream test4(resultfile.c_str());
  ofstream shb08file;
  if(test4.good()) {
    shb08file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb08file.open(resultfile.c_str());
  }
  shb08file << now << "\t" << shb08 << endl;
  shb08file.close();

  
  resultfile = "Monitoring/SHB09LeakageCurrent.dat";
  ifstream test5(resultfile.c_str());
  ofstream shb09file;
  if(test5.good()) {
    shb09file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb09file.open(resultfile.c_str());
  }
  shb09file << now << "\t" << shb09 << endl;
  shb09file.close();

  
  resultfile = "Monitoring/SHB10LeakageCurrent.dat";
  ifstream test6(resultfile.c_str());
  ofstream shb10file;
  if(test6.good()) {
    shb10file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb10file.open(resultfile.c_str());
  }
  shb10file << now << "\t" << shb10 << endl;
  shb10file.close();

  resultfile = "Monitoring/SHB11LeakageCurrent.dat";
  ifstream test0(resultfile.c_str());
  ofstream shb11file;
  if(test0.good()) {
    shb11file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb11file.open(resultfile.c_str());
  }
  shb11file << now << "\t" << shb11 << endl;
  shb11file.close();

  resultfile = "Monitoring/SHB12LeakageCurrent.dat";
  ifstream test7(resultfile.c_str());
  ofstream shb12file;
  if(test7.good()) {
    shb12file.open(resultfile.c_str(),ios::app);
  }
  else {
    shb12file.open(resultfile.c_str());
  }
  shb12file << now << "\t" << shb12 << endl;
  shb12file.close();

  resultfile = "Monitoring/SHQ13LeakageCurrent.dat";
  ifstream test8(resultfile.c_str());
  ofstream shq13file;
  if(test8.good()) {
    shq13file.open(resultfile.c_str(),ios::app);
  }
  else {
    shq13file.open(resultfile.c_str());
  }
  shq13file << now << "\t" << shq13 << endl;
  shq13file.close();

  resultfile = "Monitoring/SHQ14LeakageCurrent.dat";
  ifstream test9(resultfile.c_str());
  ofstream shq14file;
  if(test9.good()) {
    shq14file.open(resultfile.c_str(),ios::app);
  }
  else {
    shq14file.open(resultfile.c_str());
  }
  shq14file << now << "\t" << shq14 << endl;
  shq14file.close();
  
  resultfile = "Monitoring/SHQ15LeakageCurrent.dat";
  ifstream test10(resultfile.c_str());
  ofstream shq15file;
  if(test10.good()) {
    shq15file.open(resultfile.c_str(),ios::app);
  }
  else {
    shq15file.open(resultfile.c_str());
  }
  shq15file << now << "\t" << shq15 << endl;
  shq15file.close();

  resultfile = "Monitoring/SHQ16LeakageCurrent.dat";
  ifstream test11(resultfile.c_str());
  ofstream shq16file;
  if(test11.good()) {
    shq16file.open(resultfile.c_str(),ios::app);
  }
  else {
    shq16file.open(resultfile.c_str());
  }
  shq16file << now << "\t" << shq16 << endl;
  shq16file.close();
  
  resultfile = "Monitoring/pad1LeakageCurrent.dat";
  ifstream test12(resultfile.c_str());
  ofstream pad1file;
  if(test12.good()) {
    pad1file.open(resultfile.c_str(),ios::app);
  }
  else {
    pad1file.open(resultfile.c_str());
  }
  pad1file << now << "\t" << pad1 << endl;
  pad1file.close();
  
  resultfile = "Monitoring/pad2LeakageCurrent.dat";
  ifstream test13(resultfile.c_str());
  ofstream pad2file;
  if(test13.good()) {
    pad2file.open(resultfile.c_str(),ios::app);
  }
  else {
    pad2file.open(resultfile.c_str());
  }
  pad2file << now << "\t" << pad2 << endl;
  pad2file.close();
  
  resultfile = "Monitoring/pad3LeakageCurrent.dat";
  ifstream test14(resultfile.c_str());
  ofstream pad3file;
  if(test14.good()) {
    pad3file.open(resultfile.c_str(),ios::app);
  }
  else {
    pad3file.open(resultfile.c_str());
  }
  pad3file << now << "\t" << pad3 << endl;
  pad3file.close();

  resultfile = "Monitoring/pad4LeakageCurrent.dat";
  ifstream test15(resultfile.c_str());
  ofstream pad4file;
  if(test15.good()) {
    pad4file.open(resultfile.c_str(),ios::app);
  }
  else {
    pad4file.open(resultfile.c_str());
  }
  pad4file << now << "\t" << pad4 << endl;
  pad4file.close();

  resultfile = "Monitoring/TrificPositiveLeakageCurrent.dat";
  ifstream test16(resultfile.c_str());
  ofstream trifPfile;
  if(test16.good()) {
    trifPfile.open(resultfile.c_str(),ios::app);
  }
  else {
    trifPfile.open(resultfile.c_str());
  }
  trifPfile << now << "\t" << trifPos << endl;
  trifPfile.close();

  resultfile = "Monitoring/TrificNegativeLeakageCurrent.dat";
  ifstream test17(resultfile.c_str());
  ofstream trifNfile;
  if(test17.good()) {
    trifNfile.open(resultfile.c_str(),ios::app);
  }
  else {
    trifNfile.open(resultfile.c_str());
  }
  trifNfile << now << "\t" << trifNeg << endl;
  trifNfile.close();

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(4,5);
  c1->cd();

  c1->cd(1);
  TGraph *gr = new TGraph("Monitoring/SHB05LeakageCurrent.dat");
  gr->Draw("A*");
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr->SetTitle("SHB05 Leakage Current");

  c1->cd(2);
  TGraph *gr2 = new TGraph("Monitoring/SHB06LeakageCurrent.dat");
  gr2->Draw("A*");
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr2->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr2->SetTitle("SHB06 Leakage Current");

  c1->cd(3);
  TGraph *gr3 = new TGraph("Monitoring/SHB07LeakageCurrent.dat");
  gr3->Draw("A*");
  gr3->GetXaxis()->SetTimeDisplay(1);
  gr3->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr3->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr3->SetTitle("SHB07 Leakage Current");
  
  c1->cd(4);
  TGraph *gr4 = new TGraph("Monitoring/SHB08LeakageCurrent.dat");
  gr4->Draw("A*");
  gr4->GetXaxis()->SetTimeDisplay(1);
  gr4->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr4->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr4->SetTitle("SHB08 Leakage Current");

  c1->cd(5);
  TGraph *gr5 = new TGraph("Monitoring/SHB09LeakageCurrent.dat");
  gr5->Draw("A*");
  gr5->GetXaxis()->SetTimeDisplay(1);
  gr5->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr5->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr5->SetTitle("SHB09 Leakage Current");

  c1->cd(6);
  TGraph *gr6 = new TGraph("Monitoring/SHB10LeakageCurrent.dat");
  gr6->Draw("A*");
  gr6->GetXaxis()->SetTimeDisplay(1);
  gr6->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr6->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr6->SetTitle("SHB10 Leakage Current");

  c1->cd(7);
  TGraph *gr7 = new TGraph("Monitoring/SHB11LeakageCurrent.dat");
  gr7->Draw("A*");
  gr7->GetXaxis()->SetTimeDisplay(1);
  gr7->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr7->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr7->SetTitle("SHB11 Leakage Current");

  c1->cd(8);
  TGraph *gr8 = new TGraph("Monitoring/SHB12LeakageCurrent.dat");
  gr8->Draw("A*");
  gr8->GetXaxis()->SetTimeDisplay(1);
  gr8->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr8->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr8->SetTitle("SHB12 Leakage Current");

  c1->cd(9);
  TGraph *gr9 = new TGraph("Monitoring/SHQ13LeakageCurrent.dat");
  gr9->Draw("A*");
  gr9->GetXaxis()->SetTimeDisplay(1);
  gr9->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr9->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr9->SetTitle("SHQ13 Leakage Current");

  c1->cd(10);
  TGraph *gr10 = new TGraph("Monitoring/SHQ14LeakageCurrent.dat");
  gr10->Draw("A*");
  gr10->GetXaxis()->SetTimeDisplay(1);
  gr10->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr10->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr10->SetTitle("SHQ14 Leakage Current");

  c1->cd(11);
  TGraph *gr11 = new TGraph("Monitoring/SHQ15LeakageCurrent.dat");
  gr11->Draw("A*");
  gr11->GetXaxis()->SetTimeDisplay(1);
  gr11->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr11->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr11->SetTitle("SHQ15 Leakage Current");

  c1->cd(12);
  TGraph *gr12 = new TGraph("Monitoring/SHQ16LeakageCurrent.dat");
  gr12->Draw("A*");
  gr12->GetXaxis()->SetTimeDisplay(1);
  gr12->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr12->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr12->SetTitle("SHQ16 Leakage Current");


  c1->cd(13);
  TGraph *gr13 = new TGraph("Monitoring/pad1LeakageCurrent.dat");
  gr13->Draw("A*");
  gr13->GetXaxis()->SetTimeDisplay(1);
  gr13->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr13->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr13->SetTitle("Pad 1 Leakage Current");

  c1->cd(14);
  TGraph *gr14 = new TGraph("Monitoring/pad2LeakageCurrent.dat");
  gr14->Draw("A*");
  gr14->GetXaxis()->SetTimeDisplay(1);
  gr14->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr14->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr14->SetTitle("Pad 2 Leakage Current");

  c1->cd(15);
  TGraph *gr15 = new TGraph("Monitoring/pad3LeakageCurrent.dat");
  gr15->Draw("A*");
  gr15->GetXaxis()->SetTimeDisplay(1);
  gr15->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr15->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr15->SetTitle("Pad 3 Leakage Current");

  c1->cd(16);
  TGraph *gr16 = new TGraph("Monitoring/pad4LeakageCurrent.dat");
  gr16->Draw("A*");
  gr16->GetXaxis()->SetTimeDisplay(1);
  gr16->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr16->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr16->SetTitle("Pad 4 Leakage Current");


  c1->cd(17);
  TGraph *gr17 = new TGraph("Monitoring/TrificPositiveLeakageCurrent.dat");
  gr17->Draw("A*");
  gr17->GetXaxis()->SetTimeDisplay(1);
  gr17->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr17->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr17->SetTitle("TRIFIC Positive Leakage Current");

  c1->cd(18);
  TGraph *gr18 = new TGraph("Monitoring/TrificNegativeLeakageCurrent.dat");
  gr18->Draw("A*");
  gr18->GetXaxis()->SetTimeDisplay(1);
  gr18->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr18->GetYaxis()->SetTitle("Leakage Current (#muA)");
  gr18->SetTitle("TRIFIC Negative Leakage Current");


  c1->Draw();
  c1->Modified();
  c1->Update();
  c1->Draw();


 /* TCanvas *c1 = new TCanvas();
  c1->Draw();
  c1->cd();
  TGraph *gr = new TGraph("Monitoring/SHB05LeakageCurrent.dat");
  //gr->Draw("A*");
  gr->GetXaxis()->SetTimeDisplay(1);
  gr->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr->GetYaxis()->SetTitle("Leakage Current (#muA)");

  TCanvas *c2 = new TCanvas();
  c2->Draw();
  c2->cd();
  TGraph *gr2 = new TGraph("Monitoring/SHB06LeakageCurrent.dat");
  gr2->Draw("A*");
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr2->GetYaxis()->SetTitle("Leakage Current (#muA)");

  TCanvas *c3 = new TCanvas();
  c->Draw();
  c2->cd();
  TGraph *gr2 = new TGraph("Monitoring/SHB06LeakageCurrent.dat");
  gr2->Draw("A*");
  gr2->GetXaxis()->SetTimeDisplay(1);
  gr2->GetXaxis()->SetTimeFormat("#splitline{%y-%m-%d}{%H:%M:%S}%F1970-01-01 00:00:00");
  gr2->GetYaxis()->SetTitle("Leakage Current (#muA)");




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
*/
}
