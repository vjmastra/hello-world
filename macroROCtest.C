#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <TCanvas.h>
#include <TVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TH1.h>
#include <TPaveText.h>
#define NMAX 50

using namespace std;

int macroROCtest(){

  TString path1 = "./dalitzROC.txt";
  TString path2 = "./angROC.txt";
 
  ifstream myFile1(path1.Data());
  ifstream myFile2(path2.Data());

  if(!myFile1.is_open()) {
    cout << "File1 not found. Returning." << endl;
    return 0;
  }
  if(!myFile2.is_open()) {
    cout << "File2 not found. Returning." << endl;
    return 0;
  }

  vector<double> thres1, nExBins1, gof1;
  vector<double> thres2, nExBins2, gof2;
  double x0, x1, x2;

  int n1 = 0;
  do {
  myFile1 >> x0 >> x1 >> x2;
  thres1.push_back(x0);
  gof1.push_back(x1);
  nExBins1.push_back(x2);
  n1++;  
  } while (myFile1.good());

  int n2 = 0;
  do {
  myFile2 >> x0 >> x1 >> x2;
  thres2.push_back(x0);
  gof2.push_back(x1);
  nExBins2.push_back(x2);
  n2++;
  } while (myFile2.good());

  myFile1.close();
  myFile2.close();

  TGraph* massRocGraph = new TGraph(n1, nExbins1, gof1);
  TGraph* angRocGraph = new TGraph(nRes, nExbins2, gof2);
  massRocGraph->SetMarkerStyle(7);
  angRocGraph->SetMarkerStyle(7);

  TLine* lineHor = new TLine(0, +0.9, 300, +0.9);

  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("c0", "gof Dalitz Efficiency", 800, 600);
  massRocGraph->Draw("AP");
  lineHor->Draw("same");

  TCanvas *c1 = new TCanvas("c1", "gof Angle Efficiency", 800, 600);
  angRocGraph->Draw("AP");
  lineHor->Draw("same");

  return 0;
}
