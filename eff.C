#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <TString.h>
#include <TFile.h>
#include <TH2F.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMath.h>

#define NORD 3 
#define NORDANG 4 

using namespace std;

double cheb(double x, int ord){

  double val;

  switch (ord) {
  case 0: val = 1;
    break;
  case 1: val = x;
    break;
  case 2: val = 2*pow(x, 2) - 1;
    break;
  case 3: val = 4*pow(x, 3) - 3*x;
    break;
  case 4: val = 8*pow(x, 4) - 8*pow(x, 2) + 1;
    break;
  case 5: val = 16*pow(x, 5) - 20*pow(x, 3) + 5*x;
    break;
  case 6: val = 32*pow(x, 6) - 48*pow(x, 4) + 18*pow(x, 2) - 1;
    break;
  case 7: val = 64*pow(x, 7) - 112*pow(x, 5) + 56*pow(x, 3) - 7*x;
    break;
  case 8: val = 128*pow(x, 8) - 256*pow(x, 6) + 160*pow(x, 4) - 32*pow(x, 2) + 1;
    break;
  case 9: val = 256*pow(x, 9) - 576*pow(x, 7) + 432*pow(x, 5) - 120*pow(x, 3) + 9*x;
    break;
  default: val = 2*x*cheb(x, ord-1) - cheb(x, ord-2); 
  }

  return val;
}

double ellEff(double *x, double* par){

  double val = 0;
  double val1, val2;

  val1 = pow((x[0] - par[0])/par[2], 2) + pow((x[1] - par[1])/par[3], 2) - 1;
  val1 = par[4]*val1;
  val1 = 0.5*(1-erf(val1));

  int n1 = 5;
  val2 = 0;
  for (int i = 0; i < NORD+1; i++) {
    for (int j = 0; j < NORD+1; j++) {
      val2 += par[n1 + i + j*(NORD+1)]*pow(x[0], i)*pow(x[1], j);
    }
  }
  val = val1*val2;
  return val;
}

double dalitzEff(double *x, double* par){

  double MKaon = 0.493677; double MKaon2 = MKaon*MKaon;
  double MPion = 0.13957018; double MPion2 = MPion*MPion;

  double MBd = 5.27961; double MBd2 = MBd*MBd;
  double MPsi_nS = 3.096916; double MJpsi2 = MPsi_nS*MPsi_nS;

  double mKP_1 = x[0];
  double mPsiP_1 = x[1];

  double val1, val2;
  if ((mKP_1 < MKaon + MPion) || (mKP_1 > MBd - MPsi_nS) || (mPsiP_1 < MPsi_nS + MPion) || (mPsiP_1 > MBd - MKaon))
    val1 = 0;
  else { // Dalitz border from PDG KINEMATICS 43.4.3.1.
    Float_t E_P = (mPsiP_1*mPsiP_1 - MJpsi2 + MPion2)/(2*mPsiP_1) ;
    Float_t E_K = (MBd2 - mPsiP_1*mPsiP_1 - MKaon2)/(2*mPsiP_1) ;
    Float_t E_PpE_K_2 = TMath::Power((E_P + E_K),2);
    Float_t sqrt_E_P2mMP2 = TMath::Sqrt(E_P*E_P - MPion2);
    Float_t sqrt_E_K2mMK2 = TMath::Sqrt(E_K*E_K - MKaon2);
    Float_t mKP2_min = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 + sqrt_E_K2mMK2,2);
    Float_t mKP2_max = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 - sqrt_E_K2mMK2,2);
    if ((mKP_1*mKP_1 < mKP2_min) || (mKP_1*mKP_1 > mKP2_max))
      val1 = 0;
    else val1 = 1;
  }
 
  val2 = 0;
  /*
    for (int i = 0; i < NORD+1; i++) {
    for (int j = 0; j < NORD+1; j++) {
    val2 += par[i + j*(NORD+1)]*pow(x[0], i)*pow(x[1], j);
    }
    }
  */
  for (int i = 0; i < NORD+1; i++) {
    for (int j = 0; j < NORD+1; j++) {
      val2 += par[i + j*(NORD+1)]*cheb(x[0], i)*cheb(x[1], j);
    }
  }

  val = val1*val2;
  return val;
}

double rectEff(double *x, double *par){

  double val = 0;
  double valx, valy, val1, val2;
  double xmin, xmax, ymin, ymax;

  valx = 1; valy = 1;
  xmin = -1; xmax = -xmin;
  ymin = -TMath::Pi(); ymax = -ymin;

  //  valx = 0.25*(1 + erf(par[0]*(x[0] - xmin)))*(1 - erf(par[0]*(x[0] - xmax)));
  //  valy = 0.25*(1 + erf(par[1]*(x[1] - ymin)))*(1 - erf(par[1]*(x[1] - ymax)));
  if (x[0] < xmin || x[0] > xmax) valx = 0;
  if (x[1] < ymin || x[1] > ymax) valy = 0;
  val1 = valx*valy;

  int n1 = 2;
  val2 = 0;
  /*
    for (int i = 0; i < NORDANG+1; i++) {
    for (int j = 0; j < NORDANG+1; j++) {
    val2 += par[n1 + i + j*(NORDANG+1)]*pow(x[0], i)*pow(x[1], j);
    }
    }
  */

  for (int i = 0; i < NORDANG+1; i++) {
    for (int j = 0; j < NORDANG+1; j++) {
      val2 += par[n1 + i + j*(NORDANG+1)]*cheb(x[0], i)*cheb(x[1], j);
    }
  }

  val = val1*val2;
  return val;
}

int eff(int ellEffMass = 0){

  gStyle->SetOptStat(00000);

  TFile* effFile = 0;
  TString path = "/lustrehome/cristella/work/Z_analysis/exclusive/clean_14ott/original/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/";

  TString effName = path + "officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_6p0to9p0SB.root";

  effFile = TFile::Open(effName);
  if (!effFile) {
    cout <<"ERROR! Unable to open efficiency file \"" <<effName <<"\".\nReturning" <<endl;
    return -1;
  } else
    cout <<"\nUsing \"" <<effName <<"\" to compute efficiency correction" <<endl;

  TString relEffNameMass = "RelEff_psi2SPi_vs_KPi_B0constr_1B0";
  TString relEffNameAng = "RelEff_planesAngle_vs_cos_psi2S_helicityAngle";
  relEffNameMass.ReplaceAll("B0constr_1B0","hardCuts_1B0"); relEffNameAng.Append("_hardCuts_1B0");

  TH2F* relEffTH2Mass = (TH2F*)effFile->Get(relEffNameMass) ;
  TH2F* relEffTH2Ang = (TH2F*)effFile->Get(relEffNameAng) ;
   
  if (!(relEffTH2Mass)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameMass <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
  }
  if (!(relEffTH2Ang)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameAng <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
  }

  double xMinMass = relEffTH2Mass->GetXaxis()->GetBinLowEdge(1);
  double xMaxMass = relEffTH2Mass->GetXaxis()->GetBinUpEdge(relEffTH2Mass->GetNbinsX());
  double yMinMass = relEffTH2Mass->GetYaxis()->GetBinLowEdge(1);
  double yMaxMass = relEffTH2Mass->GetYaxis()->GetBinUpEdge(relEffTH2Mass->GetNbinsY());
  double xMinAng = relEffTH2Ang->GetXaxis()->GetBinLowEdge(1);
  double xMaxAng = relEffTH2Ang->GetXaxis()->GetBinUpEdge(relEffTH2Ang->GetNbinsX());
  double yMinAng = relEffTH2Ang->GetYaxis()->GetBinLowEdge(1);
  double yMaxAng = relEffTH2Ang->GetYaxis()->GetBinUpEdge(relEffTH2Ang->GetNbinsY());

  //Mass efficiency

  int n1 = 0;
  const int n2 = pow(NORD+1, 2);

  if (ellEffMass) n1 = 5;

  double par[n1+n2];

  if (ellEffMass) {
    double ellPar[5] = {1.4, 4, 1, 1, 6};
    for (int index = 0; index < n1; index++)
      par[index] = ellPar[index];
  } 

  for (int index = n1; index < n1+n2; index++)
    par[index] = 0.0;

  par[0] = 5; //c00 = 1
  par[1] = -3;
  par[2] = 2;
  par[3] = -0.1;
  par[4] = -1;
  par[5] = -1;

  TF2* fitFun;

  if (ellEffMass)
    fitFun = new TF2("fitMassEff", ellEff, xMinMass, xMaxMass, yMinMass, yMaxMass, n1+n2);
  else fitFun = new TF2("fitMassEff", dalitzEff, xMinMass, xMaxMass, yMinMass, yMaxMass, n1+n2);

  fitFun->SetParameters(par);

  double margin = 5;

  for (int index = 0; index < n1+n2; index++)
    fitFun->SetParLimits(index, par[index]-margin, par[index]+margin);

  relEffTH2Mass->Fit("fitMassEff");
 
  const int nbins = 200;

  TH1F* xProj = new TH1F("mKPiEff", "mKPiEff", relEffTH2Mass->GetNbinsX(), xMinMass, xMaxMass);
  TH1F* yProj = new TH1F("mPsiPiEff", "mPsiPiEff", relEffTH2Mass->GetNbinsY(), yMinMass, yMaxMass);
  TH1F* xFitProj = new TH1F("mKPiEffFit","mKPiEffFit", nbins, xMinMass, xMaxMass);
  TH1F* yFitProj = new TH1F("mPsiPiEffFit","mPsiPiEffFit", nbins, yMinMass, yMaxMass);

  double sum = 0;

  for (int i = 0; i < xProj->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffTH2Mass->GetNbinsY(); j++) {
      sum += relEffTH2Mass->GetBinContent(i+1, j+1);
    }
    sum *= (yMaxMass - yMinMass)/relEffTH2Mass->GetNbinsY();
    xProj->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < yProj->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < relEffTH2Mass->GetNbinsX(); i++) {
      sum += relEffTH2Mass->GetBinContent(i+1, j+1);
    }
    sum *= (xMaxMass - xMinMass)/relEffTH2Mass->GetNbinsX();
    yProj->SetBinContent(j+1, sum);
  }

  double deltax = (xMaxMass - xMinMass)/nbins;
  double deltay = (yMaxMass - yMinMass)/nbins;
  double xvar, yvar;

  for (int i = 0; i < nbins; i++) {
    sum = 0;
    xvar = xMinMass + (i + 0.5)*deltax;
    for (int j = 0; j < nbins; j++) {
      yvar = yMinMass + (j + 0.5)*deltay;
      sum += fitFun->Eval(xvar, yvar);
    }
    sum *= deltay;
    xFitProj->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < nbins; j++) {
    sum = 0;
    yvar = yMinMass + (j + 0.5)*deltay;
    for (int i = 0; i < nbins; i++) {
      xvar = xMinMass + (i + 0.5)*deltax;
      sum += fitFun->Eval(xvar, yvar);
    }
    sum *= deltax;
    yFitProj->SetBinContent(j+1, sum);
  }

  //Angle efficiency

  const int n1a = 0; ///perchè ho tolto smoothing erf, con erf = 2
  const int n2a = pow(NORDANG+1, 2);
  double rectPar[2] = {6, 6}; 
  double parAng[n1a+n2a]; 

  for (int index = 0; index < n1a; index++)
    parAng[index] = rectPar[index];

  parAng[n1a] = 1; //c00 = 1
  for (int index = n1a+1; index < n1a+n2a; index++)
    parAng[index] = 0;

  TF2* fitFunAng;
  fitFunAng = new TF2("fitAngEff", rectEff, xMinAng, xMaxAng, yMinAng, yMaxAng, n1a+n2a);
  fitFunAng->SetParameters(parAng);

  for (int index = 0; index < n1a+n2a; index++)
    fitFun->SetParLimits(index, parAng[index]-margin, parAng[index]+margin);

  relEffTH2Ang->Fit("fitAngEff");

  TH1F* xProjAng = new TH1F("cosMuMuEff", "cosMuMuEff", relEffTH2Ang->GetNbinsX(), xMinAng, xMaxAng);
  TH1F* yProjAng = new TH1F("phiEff", "phiEff", relEffTH2Ang->GetNbinsY(), yMinAng, yMaxAng);
  TH1F* xFitProjAng = new TH1F("cosMuMuEffFit","cosMuMuEffFit", nbins, xMinAng, xMaxAng);
  TH1F* yFitProjAng = new TH1F("phiEffFit","phiEffFit", nbins, yMinAng, yMaxAng);

  sum = 0;

  for (int i = 0; i < xProjAng->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffTH2Ang->GetNbinsY(); j++) {
      sum += relEffTH2Ang->GetBinContent(i+1, j+1);
    }
    sum *= (yMaxAng - yMinAng)/relEffTH2Ang->GetNbinsY();
    xProjAng->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < yProjAng->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < relEffTH2Ang->GetNbinsX(); i++) {
      sum += relEffTH2Ang->GetBinContent(i+1, j+1);
    }
    sum *= (xMaxAng - xMinAng)/relEffTH2Ang->GetNbinsX();
    yProjAng->SetBinContent(j+1, sum);
  }

  deltax = (xMaxAng - xMinAng)/nbins;
  deltay = (yMaxAng - yMinAng)/nbins;

  for (int i = 0; i < nbins; i++) {
    sum = 0;
    xvar = xMinAng + (i + 0.5)*deltax;
    for (int j = 0; j < nbins; j++) {
      yvar = yMinAng + (j + 0.5)*deltay;
      sum += fitFunAng->Eval(xvar, yvar);
    }
    sum *= deltay;
    xFitProjAng->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < nbins; j++) {
    sum = 0;
    yvar = yMinAng + (j + 0.5)*deltay;
    for (int i = 0; i < nbins; i++) {
      xvar = xMinAng + (i + 0.5)*deltax;
      sum += fitFunAng->Eval(xvar, yvar);
    }
    sum *= deltax;
    yFitProjAng->SetBinContent(j+1, sum);
  }

  //Plot

  relEffTH2Mass->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  relEffTH2Mass->GetYaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  relEffTH2Ang->GetXaxis()->SetTitle("cos(MuMu)");
  relEffTH2Ang->GetYaxis()->SetTitle("phi [rad]");

  xProj->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  yProj->GetXaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  xProjAng->GetXaxis()->SetTitle("cos(MuMu)");
  yProjAng->GetXaxis()->SetTitle("phi [rad]");

  xFitProj->SetLineStyle(2);
  xFitProj->SetLineColor(2);
  xFitProj->SetLineWidth(1);  
  yFitProj->SetLineStyle(2);
  yFitProj->SetLineColor(2);
  yFitProj->SetLineWidth(1);
  xFitProjAng->SetLineStyle(2);
  xFitProjAng->SetLineColor(2);
  xFitProjAng->SetLineWidth(1);
  yFitProjAng->SetLineStyle(2);
  yFitProjAng->SetLineColor(2);
  yFitProjAng->SetLineWidth(1);
/*
  TCanvas *c0 = new TCanvas("canvas", "canvas", 800, 600);
  c0->cd();
  relEffTH2Mass->Draw("COLZ");
  c0->SaveAs("massEff.png");

  TCanvas *c1 = new TCanvas("canvas1", "canvas1", 800, 600);
  c1->cd();
  relEffTH2Ang->Draw("COLZ");
  c1->SaveAs("angEff.png");
*/ 
  TCanvas *cxMass = new TCanvas("canvasDpProjX", "canvasDpProjX", 800, 600);
  cxMass->cd();
  xProj->Draw();
  xFitProj->Draw("same");
  xProj->Draw("same");

  TCanvas *cyMass = new TCanvas("canvasDpProjY", "canvasDpProjY", 800, 600);  
  cyMass->cd();
  yProj->Draw();
  yFitProj->Draw("same");
  yProj->Draw("same");

  TCanvas *cxAng = new TCanvas("canvasAngProjX", "canvasAngProjX", 800, 600);
  cxAng->cd();
  xProjAng->Draw();
  xFitProjAng->Draw("same");
  xProjAng->Draw("same");

  TCanvas *cyAng = new TCanvas("canvasAngProjY", "canvasAngProjY", 800, 600);
  cyAng->cd();
  yProjAng->Draw();
  yFitProjAng->Draw("same");
  yProjAng->Draw("same");

  return 0;

}
