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

  double thres1[NMAX], nExBins1[NMAX], gof1[NMAX];
  double thres2[NMAX], nExBins2[NMAX], gof2[NMAX];
  double x0, x1, x2;

  int n1 = 0;
  for (int i = 0; i < NMAX; i++) {
  if (myFile1.good())
     myFile1 >> x0 >> x1 >> x2;
  else {
    x0 = -1; x1 = -1; x2 = 0;
    }
  thres1[i] = x0;
  gof1[i] = x1;
  nExBins1[i] = x2;
  n1++;  
  }

  int n2 = 0;
  for (int i = 0; i < NMAX; i++) {
  if (myFile2.good())
  myFile2 >> x0 >> x1 >> x2;
  else {
    x0 = -1; x1 = -1; x2 = 0;
    }
  thres2[i] = x0;
  gof2[i] = x1;
  nExBins2[i] = x2;
  n2++;
  } while (myFile2.good());

  myFile1.close();
  myFile2.close();

  double gofVal = 0.85;
  double inf, sup;
  int i1inf, i1sup;
  int i2inf, i2sup;

  inf = -0.1, sup = 1.1;
  for (int i = 0; i < NMAX; i++) {
  if (gof1[i] < gofVal && gof1[i] > inf) {
  i1inf = i;
  inf = gof1[i];
  } 
  if (gof1[i] > gofVal && gof1[i] < sup) {
  i1sup = i;
  sup = gof1[i];
  }
  } 
  
  inf = -0.1, sup = 1.1;
  for (int i = 0; i < NMAX; i++) {
  if (gof2[i] < gofVal && gof2[i] > inf) {
  i2inf = i;
  inf = gof2[i];
  } 
  if (gof2[i] > gofVal && gof2[i] < sup) {
  i2sup = i;
  sup = gof2[i];
  }
  } 

  double thMass, thAng;
  double nMass, nAng;

  thMass = 0.5*(thres1[i1inf] + thres1[i1sup]);
  thAng = 0.5*(thres2[i2inf] + thres2[i2sup]);
  nMass = 0.5*(nExBins1[i1inf] + nExBins1[i1sup]);
  nAng = 0.5*(nExBins2[i2inf] + nExBins2[i2sup]);

  cout << "Approx. threshold on Chi2 for a (1-p) = " << gofVal << endl;
  cout << "Dalitz: " << thMass << " excluding " << nMass << " bins" << endl;
  cout << "Angles: " << thAng << " excluding " << nAng << " bins" << endl;

  TGraph* massRocGraph = new TGraph(NMAX, nExBins1, gof1);
  TGraph* angRocGraph = new TGraph(NMAX, nExBins2, gof2);
  massRocGraph->SetMarkerStyle(7);
  massRocGraph->SetTitle("Dalitz Efficiency ROC curve");
  massRocGraph->GetXaxis()->SetTitle("#bins excluded");
  massRocGraph->GetYaxis()->SetTitle("g.o.f");
  massRocGraph->GetXaxis()->SetRangeUser(0, 250);
  massRocGraph->GetYaxis()->SetRangeUser(0, 1.1);
  angRocGraph->SetMarkerStyle(7);
  angRocGraph->SetTitle("Angles Efficiency ROC curve");
  angRocGraph->GetXaxis()->SetTitle("#bins excluded");
  angRocGraph->GetYaxis()->SetTitle("g.o.f");
  angRocGraph->GetXaxis()->SetRangeUser(0, 250);
  angRocGraph->GetYaxis()->SetRangeUser(0, 1.1);

  TLine* lineHor = new TLine(0, gofVal, 210, gofVal);
  lineHor->SetLineStyle(2);
  lineHor->SetLineColor(2);
  TLine* lineHor2 = new TLine(0, gofVal, 235, gofVal);
  lineHor2->SetLineStyle(2);
  lineHor2->SetLineColor(2);

  gStyle->SetOptStat(0);
  TCanvas *c0 = new TCanvas("c0", "dalitzRoc", 800, 600);
  massRocGraph->Draw("AP");
  lineHor2->Draw("same");

  TCanvas *c1 = new TCanvas("c1", "anglesRoc", 800, 600);
  angRocGraph->Draw("AP");
  lineHor->Draw("same");

  return 0;
}
