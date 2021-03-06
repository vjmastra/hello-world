#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF2.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMath.h>

#define NORD0 4
#define NORD1 4
#define NORD2 5
#define NORD3 5

using namespace std;

//cose da fare
//fill del taglio
//lettura dei valori iniziali dal file (con flag per read o readandwrite)
//flag per proiezioni con/senza scarti
//controllo su flag
//count su pdf delle proiezioni
//pulire mettendo le masse come variabili "globali"
//ripetere curva ROc con polinomi di ordine stabile
//grafica degli hist

double inDalitz(double x, double y){

  double MKaon = 0.493677; double MKaon2 = MKaon*MKaon;
  double MPion = 0.13957018; double MPion2 = MPion*MPion;

  double MBd = 5.27961; double MBd2 = MBd*MBd;
  double MPsi_nS = 3.096916; double MJpsi2 = MPsi_nS*MPsi_nS;

  double mKP_1 = x;
  double mPsiP_1 = y;

  if ((mKP_1 < MKaon + MPion) || (mKP_1 > MBd - MPsi_nS) || (mPsiP_1 < MPsi_nS + MPion) || (mPsiP_1 > MBd - MKaon))
    return 0;
  else { // Dalitz border from PDG KINEMATICS 43.4.3.1.
    Float_t E_P = (mPsiP_1*mPsiP_1 - MJpsi2 + MPion2)/(2*mPsiP_1) ;
    Float_t E_K = (MBd2 - mPsiP_1*mPsiP_1 - MKaon2)/(2*mPsiP_1) ;
    Float_t E_PpE_K_2 = TMath::Power((E_P + E_K),2);
    Float_t sqrt_E_P2mMP2 = TMath::Sqrt(E_P*E_P - MPion2);
    Float_t sqrt_E_K2mMK2 = TMath::Sqrt(E_K*E_K - MKaon2);
    Float_t mKP2_min = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 + sqrt_E_K2mMK2,2);
    Float_t mKP2_max = E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 - sqrt_E_K2mMK2,2);
    if ((mKP_1*mKP_1 < mKP2_min) || (mKP_1*mKP_1 > mKP2_max))
      return 0;
    else return 1;
  }

}

double xDalitz(double y, int side){

  //There are two x-values for each y-value
  //The x-values are labeled: 0-left, 1-right

  if (side != 0 && side != 1) return -1;

  double MKaon = 0.493677; double MKaon2 = MKaon*MKaon;
  double MPion = 0.13957018; double MPion2 = MPion*MPion;

  double MBd = 5.27961; double MBd2 = MBd*MBd;
  double MPsi_nS = 3.096916; double MJpsi2 = MPsi_nS*MPsi_nS;

  double mPsiP_1 = y;

  if ((mPsiP_1 < MPsi_nS + MPion) || (mPsiP_1 > MBd - MKaon))
    return -1;
  else { // Dalitz border from PDG KINEMATICS 43.4.3.1.
    Float_t E_P = (mPsiP_1*mPsiP_1 - MJpsi2 + MPion2)/(2*mPsiP_1) ;
    Float_t E_K = (MBd2 - mPsiP_1*mPsiP_1 - MKaon2)/(2*mPsiP_1) ;
    Float_t E_PpE_K_2 = TMath::Power((E_P + E_K),2);
    Float_t sqrt_E_P2mMP2 = TMath::Sqrt(E_P*E_P - MPion2);
    Float_t sqrt_E_K2mMK2 = TMath::Sqrt(E_K*E_K - MKaon2);
    Float_t mKP_min = pow(E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 + sqrt_E_K2mMK2,2), 0.5);
    Float_t mKP_max = pow(E_PpE_K_2 - TMath::Power(sqrt_E_P2mMP2 - sqrt_E_K2mMK2,2), 0.5);
    if (side == 0) return mKP_min;
    else if (side == 1) return mKP_max;
  }

  return -1;

}

double cheb(double x, int ord){

  double val;
 
  if (ord < 0) val = 0;

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

int factorial(int n){

  if (n == 0) return 1;
  else return n*factorial(n-1);

}

int binomial(int k, int n){

  if (k > n) return 0;
  else return factorial(n)/factorial(k)*factorial(n-k);

}

double bernstein(double x, int i, int order, int var){

  double a, b;
  switch (var) {
  case 0: a = 0.6; b = 2.2; //MassKPi
    break;
  case 1: a = 3.2; b = 4.8; //MassPsiPi
    break;
  case 2: a = -1; b = 1; //CosMuMu
    break;
  case 3: a = -3.14; b = 3.14; //Phi
    break;
  default: return 0;

  }

  x = (x-a)/(b-a);

  return binomial(i, order)*pow(x, i)*pow(1-x, order-i);

}

double ellEff(double *x, double* par){

  double val = 0;
  double val1, val2;

  val1 = pow((x[0] - par[0])/par[2], 2) + pow((x[1] - par[1])/par[3], 2) - 1;
  val1 = par[4]*val1;
  val1 = 0.5*(1-erf(val1));

  int n1 = 5;
  val2 = 0;
  for (int i = 0; i < NORD0+1; i++) {
    for (int j = 0; j < NORD1+1; j++) {
      val2 += par[n1 + i + j*(NORD1+1)]*pow(x[0], i)*pow(x[1], j);
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
  /* x-y Polynomial
     for (int i = 0; i < NORD0+1; i++) {
     for (int j = 0; j < NORD1+1; j++) {
     val2 += par[i + j*(NORD1+1)]*pow(x[0], i)*pow(x[1], j);
     }
     }
  */
  /* Chebchev Polynomial
     for (int i = 0; i < NORD0+1; i++) {
     for (int j = 0; j < NORD1+1; j++) {
     val2 += par[i + j*(NORD1+1)]*cheb(x[0], i)*cheb(x[1], j);
     }
     }
  */
  //Bernstein Polynomial
  for (int i = 0; i < NORD0+1; i++) {
    for (int j = 0; j < NORD1+1; j++) {
      val2 += par[i + j*(NORD1+1)]*bernstein(x[0], i, NORD0, 0)*bernstein(x[1], j, NORD1, 1);
    }
  }

  return val1*(1+val2);

}

double erfDalitzEff(double *x, double *par){

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
    double xmin = TMath::Sqrt(mKP2_min);
    double xmax = TMath::Sqrt(mKP2_max);
    val1 = 0.25*(1 + erf(5*(x[0] - xmin)))*(1 - erf(5*(x[0] - xmax)));
  }

  int n1 = 0;
  val2 = 0;

  for (int i = 0; i < NORD0+1; i++) {
    for (int j = 0; j < NORD1+1; j++) {
      val2 += par[n1 + i + j*(NORD1+1)]*bernstein(x[0], i, NORD0, 0)*bernstein(x[1], j, NORD1, 1);
    }
  }

  //Chebchev Polynomial
  /*
    for (int i = 0; i < NORD0+1; i++) {
    for (int j = 0; j < NORD1+1; j++) {
    val2 += par[i + j*(NORD1+1)]*cheb(x[0], i)*cheb(x[1], j);
    }
    }
  */
  val = val1*(1+val2);
  return val;
}

double rectEff(double *x, double *par){

  double val = 0;
  double valx, valy, val1, val2;
  double xmin, xmax, ymin, ymax;

  valx = 1; valy = 1;
  xmin = -1; xmax = -xmin;
  ymin = -TMath::Pi(); ymax = -ymin;

  if (x[0] < xmin || x[0] > xmax) valx = 0;
  if (x[1] < ymin || x[1] > ymax) valy = 0;
  val1 = valx*valy;

  val2 = 0;
  /*
    for (int i = 0; i < NORD2+1; i++) {
    for (int j = 0; j < NORD3+1; j++) {
    val2 += par[i + j*(NORD3+1)]*pow(x[0], i)*pow(x[1], j);
    }
    }
  */
  /*
    for (int i = 0; i < NORD2+1; i++) {
    for (int j = 0; j < NORD3+1; j++) {
    val2 += par[i + j*(NORD3+1)]*cheb(x[0], i)*cheb(x[1], j);
    }
    }
  */
  for (int i = 0; i < NORD2+1; i++) {
    for (int j = 0; j < NORD3+1; j++) {
      val2 += par[i + j*(NORD3+1)]*bernstein(x[0], i, NORD2, 2)*bernstein(x[1], j, NORD3, 3);
    }
  }

  val = val1*(1+val2);
  return val;
}

double erfRectEff(double *x, double *par){

  double val = 0;
  double valx, valy, val1, val2;
  double xmin, xmax, ymin, ymax;

  valx = 1; valy = 1;
  xmin = -1; xmax = -xmin;
  ymin = -TMath::Pi(); ymax = -ymin;

  valx = 0.25*(1 + erf(par[0]*(x[0] - xmin)))*(1 - erf(par[0]*(x[0] - xmax)));
  valy = 0.25*(1 + erf(par[1]*(x[1] - ymin)))*(1 - erf(par[1]*(x[1] - ymax)));
  val1 = valx*valy;
  
  int n1 = 2;
  val2 = 0;

  for (int i = 0; i < NORD2+1; i++) {
    for (int j = 0; j < NORD3+1; j++) {
      val2 += par[n1 + i + j*(NORD3+1)]*bernstein(x[0], i, NORD2, 2)*bernstein(x[1], j, NORD3, 3);
    }
  }

  val = val1*(1+val2);
  return val;
}

int eff(int erfSmooth = 0, int chiThres = 1, int chiRoc = 0, int noEdge = 0, int readPar = 0){

  if (chiThres == 0 && chiRoc == 1) chiRoc = 0;
  //noEdge non attiva
  //readPar da attivare

  gStyle->SetOptStat(00000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //Reading histograms

  TFile* effFile = 0;
  //  TString path = "/lustrehome/cristella/work/Z_analysis/exclusive/clean_14ott/original/CMSSW_5_3_22/src/UserCode/MuMuPiKPAT/test/sanjay/selector/";
  TString path = "./";

  TString effName = path + "officialMC_noPtEtaCuts_JPsi_Bd2MuMuKPi_2p0Sig_6p0to9p0SB_62M.root";

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
  
  relEffNameMass.ReplaceAll("Rel","Abs"); relEffNameAng.ReplaceAll("Rel", "Abs");
  /*
    TH2F* absEffTH2Mass = (TH2F*)effFile->Get(relEffNameMass) ;
    TH2F* absEffTH2Ang = (TH2F*)effFile->Get(relEffNameAng) ;

    if (!(absEffTH2Mass)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameMass <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
    }
    if (!(absEffTH2Ang)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameAng <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
    }
  */ 
  relEffNameMass.ReplaceAll("AbsEff","RelEffErr"); relEffNameAng.ReplaceAll("AbsEff", "RelEffErr");

  TH2F* relEffErrTH2Mass = (TH2F*)effFile->Get(relEffNameMass) ;
  TH2F* relEffErrTH2Ang = (TH2F*)effFile->Get(relEffNameAng) ;

  if (!(relEffErrTH2Mass)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameMass <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
  }
  if (!(relEffErrTH2Ang)) {
    std::cout<<"ERROR! Efficiency TH2 named \'"<<relEffNameAng <<"\' not found in TFile \'" <<effFile->GetName() <<"\'.\nReturning" <<std::endl;
    return -1;
  }

  double thres = 5.2;

  double xMinMass = relEffTH2Mass->GetXaxis()->GetBinLowEdge(1);
  double xMaxMass = relEffTH2Mass->GetXaxis()->GetBinUpEdge(relEffTH2Mass->GetNbinsX());
  double yMinMass = relEffTH2Mass->GetYaxis()->GetBinLowEdge(1);
  double yMaxMass = relEffTH2Mass->GetYaxis()->GetBinUpEdge(relEffTH2Mass->GetNbinsY());
  double xMinAng = relEffTH2Ang->GetXaxis()->GetBinLowEdge(1);
  double xMaxAng = relEffTH2Ang->GetXaxis()->GetBinUpEdge(relEffTH2Ang->GetNbinsX());
  double yMinAng = relEffTH2Ang->GetYaxis()->GetBinLowEdge(1);
  double yMaxAng = relEffTH2Ang->GetYaxis()->GetBinUpEdge(relEffTH2Ang->GetNbinsY());

  TH1F* relEffErrKPi = new TH1F("errKPi", "errKPi", relEffTH2Mass->GetNbinsX(), xMinMass, xMaxMass);
  TH1F* relEffErrPsiPi = new TH1F("errPsiPi", "errPsiPi", relEffTH2Mass->GetNbinsY(), yMinMass, yMaxMass);
  TH1F* relEffErrCosMuMu = new TH1F("errCosMuMu", "errCosMuMu", relEffTH2Ang->GetNbinsX(), xMinAng, xMaxAng);
  TH1F* relEffErrPhi = new TH1F("errPhi", "errPhi", relEffTH2Ang->GetNbinsY(), yMinAng, yMaxAng); 

  TH2F* chiTH2Mass = new TH2F("chiTH2Mass", "chiTH2Mass", relEffTH2Mass->GetNbinsX(), xMinMass, xMaxMass, relEffTH2Mass->GetNbinsY(), yMinMass, yMaxMass);
  TH2F* chiTH2Ang = new TH2F("chiTH2Ang", "chiTH2Ang", relEffTH2Ang->GetNbinsX(), xMinAng, xMaxAng, relEffTH2Ang->GetNbinsY(), yMinAng, yMaxAng);

  double chiKPi[relEffTH2Mass->GetNbinsX()], chiPsiPi[relEffTH2Mass->GetNbinsY()];
  double chiCosMuMu[relEffTH2Ang->GetNbinsX()], chiPhi[relEffTH2Ang->GetNbinsY()];
  double indexKPi[relEffTH2Mass->GetNbinsX()], indexPsiPi[relEffTH2Mass->GetNbinsY()];
  double indexCosMuMu[relEffTH2Ang->GetNbinsX()], indexPhi[relEffTH2Ang->GetNbinsY()];

  double chiMass, dofMass, probMass;
  double chiAng, dofAng, probAng;
  double chi2KPi, chi2PsiPi, chi2CosMuMu, chi2Phi;
  double dofKPi, dofPsiPi, dofCosMuMu, dofPhi;
  double probKPi, probPsiPi, probCosMuMu, probPhi;

  int nExcludedMass = 0;
  int nExcludedAng = 0;

  TH2F* relEffTH2MassCut = (TH2F*)relEffTH2Mass->Clone();
  TH2F* relEffTH2AngCut = (TH2F*)relEffTH2Ang->Clone();
  TH2F* chiTH2MassCut = (TH2F*)chiTH2Mass->Clone();
  TH2F* chiTH2AngCut = (TH2F*)chiTH2Ang->Clone(); 
  TH2F* relEffErrTH2MassCut = (TH2F*)relEffErrTH2Mass->Clone();
  TH2F* relEffErrTH2AngCut = (TH2F*)relEffErrTH2Ang->Clone();

  const int npy = 100; //50k continous line
  double xDal[2*npy], yDal[2*npy];
  double step = (yMaxMass - yMinMass)/npy;

  for (int i = 0; i < npy; i++) {
    for (int j = 0; j < 2; j++) {
      yDal[2*i + j] = yMinMass + (i+0.5)*step;
      xDal[2*i + j] = xDalitz(yDal[2*i + j], j);
    }
  }
 
  TGraph* dalitzGraph = new TGraph(2*npy, xDal, yDal);
  /*
    TString dalParFile = "./dalitzPar.txt";
    TString angParFile = "./angPar.txt";
    ifstream readDalPar(darParFile.Data());
    ifstream readAngPar(angParFile.Data());

    if (readPar) {
    if(!readDalPar.is_open()) {
    cout << "File dalitzPar.txt not found. Returning." << endl;
    return 0;
    }
    if(!readAngPar.is_open()) {
    cout << "File angPar.txt not found. Returning." << endl;
    return 0;
    }
    }
  */
  TString dalRocPath = "./dalitzROC.txt";
  TString angRocPath = "./angROC.txt";
  ofstream outDalRoc;
  ofstream outAngRoc;

  if (chiRoc) {
    outDalRoc.open(dalRocPath.Data(), std::ofstream::out | std::ofstream::app);
    outAngRoc.open(angRocPath.Data(), std::ofstream::out | std::ofstream::app);

    if(!outDalRoc.is_open()) {
      cout << "File dalitzPar.txt not found. Returning." << endl;
      return 0;
    }
    if(!outAngRoc.is_open()) {
      cout << "File angPar.txt not found. Returning." << endl;
      return 0;
    }
  }

  //Mass efficiency fit

  int n1 = 0;
  const int n2 = (NORD0+1)*(NORD1+1);

  if (erfSmooth == 1) n1 = 5;

  double par[n1+n2];

  if (erfSmooth == 1) {
    double ellPar[5] = {1.4, 4, 1, 1, 6};
    for (int index = 0; index < n1; index++)
      par[index] = ellPar[index];
  } 

  //sarebbe par[n1+..]
  /*  par[0] = 5; //coeff per cheb ordine 3
      par[1] = -3;
      par[2] = 2;
      par[3] = -0.2;
      par[4] = -0.6;
      par[5] = -0.7;
  */ 
  /*  if (readPar) {
      double parVal = 0; 
      for(int i = n1; i < n1+n2; i++){
      dalParFile.good() ? dalParFile >> parVal : parVal = 0;
      par[i] = parVal;
      }
      }
      else 
  */    for (int index = n1; index < n1+n2; index++)
    par[index] = 0;

  TF2* fitFun;

  if (erfSmooth == 1) fitFun = new TF2("fitMassEff", ellEff, xMinMass, xMaxMass, yMinMass, yMaxMass, n1+n2);
  if (erfSmooth == 2) fitFun = new TF2("fitMassEff", dalitzEff, xMinMass, xMaxMass, yMinMass, yMaxMass, n1+n2);
  if (erfSmooth == 0) fitFun = new TF2("fitMassEff", erfDalitzEff, xMinMass, xMaxMass, yMinMass, yMaxMass, n1+n2);
  fitFun->SetParameters(par);

  double margin = 5;
  for (int index = 0; index < n1+n2; index++)
    fitFun->SetParLimits(index, par[index]-margin, par[index]+margin);

  relEffTH2Mass->Fit("fitMassEff");
 
  //Mass Efficiency Chi2D

  double chi = 0, sigma = 0;
  double widthx = (xMaxMass - xMinMass)/relEffTH2Mass->GetNbinsX();
  double widthy = (yMaxMass - yMinMass)/relEffTH2Mass->GetNbinsY();
  const int nbins0 = relEffTH2Mass->GetNbinsX()*10;
  const int nbins1 = relEffTH2Mass->GetNbinsY()*10;
  const int nbins2 = relEffTH2Ang->GetNbinsX()*10;
  const int nbins3 = relEffTH2Ang->GetNbinsY()*10;
  int scalex = floor(nbins0/relEffTH2Mass->GetNbinsX());
  int scaley = floor(nbins1/relEffTH2Mass->GetNbinsY());
  double deltax = (xMaxMass - xMinMass)/nbins0;
  double deltay = (yMaxMass - yMinMass)/nbins1;
  double sum = 0, val = 0;
  int count = 0;
  double xvar, yvar;

  for (int i = 0; i < chiTH2Mass->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Mass->GetNbinsY(); j++) {
      sum = 0;
      count = 0;
      for (int k = 0; k < scalex; k++) {
        for (int l = 0; l < scaley; l++) {
          xvar = xMinMass + i*widthx + (k + 0.5)*deltax;
          yvar = yMinMass + j*widthy + (l + 0.5)*deltay;
          val = fitFun->Eval(xvar, yvar);
          sum += val;
          if (val != 0) count++;
        }
      }
      count == 0 ? sum = 0 : sum = sum/count;
      chi = relEffTH2Mass->GetBinContent(i+1, j+1) - sum;
      sigma = relEffErrTH2Mass->GetBinContent(i+1, j+1);
      sigma == 0 ? chi = 0 : chi = pow(chi/sigma, 2);
      chiTH2Mass->SetBinContent(i+1, j+1, chi); 
      if (chiThres && chi > thres) {
	relEffTH2MassCut->SetBinContent(i+1, j+1, 0);
	relEffErrTH2MassCut->SetBinContent(i+1, j+1, 0);
	chi = 0;
	nExcludedMass++;
      }
      chiTH2MassCut->SetBinContent(i+1, j+1, chi);
    }
  }

  chiMass = 0;
  dofMass = 0;
  for (int i = 0; i < chiTH2Mass->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Mass->GetNbinsY(); j++) {
      val = chiTH2MassCut->GetBinContent(i+1, j+1);
      chiMass += val;
      if (val != 0) dofMass++;
    }
  }

  dofMass = dofMass - (n1+n2);
  probMass = TMath::Prob(chiMass, dofMass);

  int nBinsHist0 = 48;
  TH1F* chi2MassDist = new TH1F("chi2MassDistribution", "chi2MassDistribution", nBinsHist0, 0, 12);
  for (int i = 0; i < chiTH2Mass->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Mass->GetNbinsY(); j++) {
      chi2MassDist->Fill(chiTH2Mass->GetBinContent(i+1, j+1));
    }
  }
  chi2MassDist->AddBinContent(nBinsHist0, chi2MassDist->GetBinContent(nBinsHist0+1));

  TH1F* chi2MassDistCut = new TH1F("chi2MassDistributionWithCut", "chi2MassDistributionWithCut", nBinsHist0, 0, 12);
  for (int i = 0; i < chiTH2MassCut->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2MassCut->GetNbinsY(); j++) {
      chi2MassDistCut->Fill(chiTH2MassCut->GetBinContent(i+1, j+1));
    }
  }
  chi2MassDistCut->AddBinContent(nBinsHist0, chi2MassDistCut->GetBinContent(nBinsHist0+1));

  //Mass efficiency projections (without cuts)

  //To implement cuts: TH2Cuts in proj
  //Skip excluded bins in PDF
  //Skip excluded bins in errors

  TH1F* xProj = new TH1F("mKPiEff", "mKPiEff", relEffTH2Mass->GetNbinsX(), xMinMass, xMaxMass);
  TH1F* yProj = new TH1F("mPsiPiEff", "mPsiPiEff", relEffTH2Mass->GetNbinsY(), yMinMass, yMaxMass);
  TH1F* xFitProj = new TH1F("mKPiEffFit","mKPiEffFit", nbins0, xMinMass, xMaxMass);
  TH1F* yFitProj = new TH1F("mPsiPiEffFit","mPsiPiEffFit", nbins1, yMinMass, yMaxMass);

  for (int i = 0; i < xProj->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffTH2Mass->GetNbinsY(); j++) {
      val = relEffTH2Mass->GetBinContent(i+1, j+1);
      sum += val;
    }
    xProj->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < yProj->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < relEffTH2Mass->GetNbinsX(); i++) {
      val = relEffTH2Mass->GetBinContent(i+1, j+1);
      sum += val;
    }
    yProj->SetBinContent(j+1, sum);
  }

  for (int i = 0; i < nbins0; i++) {
    sum = 0;
    xvar = xMinMass + (i + 0.5)*deltax;
    for (int j = 0; j < nbins1; j++) {
      yvar = yMinMass + (j + 0.5)*deltay;
      val = fitFun->Eval(xvar, yvar);
      sum += val;
    }
    sum *= pow(scaley, -1);
    xFitProj->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < nbins1; j++) {
    sum = 0;
    yvar = yMinMass + (j + 0.5)*deltay;
    for (int i = 0; i < nbins0; i++) {
      xvar = xMinMass + (i + 0.5)*deltax;
      val = fitFun->Eval(xvar, yvar);
      sum += val;
    }
    sum *= pow(scalex, -1);
    yFitProj->SetBinContent(j+1, sum);
  }

  //Pull MassKPi
 
  for (int i = 0; i < relEffErrTH2Mass->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffErrTH2Mass->GetNbinsY(); j++) {
      sum += pow(relEffErrTH2Mass->GetBinContent(i+1, j+1), 2);
    }
    sum = pow(sum, 0.5);
    relEffErrKPi->SetBinContent(i+1, sum); //sostituire con vettore
  }

  chi2KPi = 0;
  dofKPi = 0;

  for (int i = 0; i < xProj->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < yProj->GetNbinsX(); j++) {
      for (int k = 0; k < scalex; k++) {
        xvar = xMinMass + i*widthx + (k + 0.5)*deltax;
        for (int l = 0; l < scaley; l++) {
          yvar = yMinMass + j*widthy + (l + 0.5)*deltay;
          sum += fitFun->Eval(xvar, yvar);
	  //applicare correzione del bordo e norm su count
        }
      }
    }
    sum *= pow(scalex*scaley, -1);
    chi = xProj->GetBinContent(i+1) - sum;
    sigma = relEffErrKPi->GetBinContent(i+1);
    sigma == 0 ? chi = 0 : chi = chi/sigma;
    indexKPi[i] = xMinMass + (i+0.5)*widthx;
    chiKPi[i] = chi;
    chi2KPi += chi*chi;
    if (chi != 0) dofKPi++; 
  }

  dofKPi = dofKPi - (NORD0+1);
  probKPi = TMath::Prob(chi2KPi, probKPi);

  //Pull MassPsiPi

  for (int j = 0; j < relEffErrTH2Mass->GetNbinsY(); j++) {
    sum = 0;
    for (int i = 0; i < relEffErrTH2Mass->GetNbinsX(); i++) {
      sum += pow(relEffErrTH2Mass->GetBinContent(i+1, j+1), 2);
    }
    sum = pow(sum, 0.5);
    relEffErrPsiPi->SetBinContent(j+1, sum);
  }

  chi2PsiPi = 0;
  dofPsiPi = 0;

  for (int j = 0; j < yProj->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < xProj->GetNbinsX(); i++) {
      for (int k = 0; k < scalex; k++) {
        xvar = xMinMass + i*widthx + (k + 0.5)*deltax;
        for (int l = 0; l < scaley; l++) {
          yvar = yMinMass + j*widthy + (l + 0.5)*deltay;
          sum += fitFun->Eval(xvar, yvar);
	  //count
        }
      }
    }
    sum *= pow(scalex*scaley, -1);
    chi = yProj->GetBinContent(j+1) - sum;
    sigma = relEffErrPsiPi->GetBinContent(j+1);
    sigma == 0 ? chi = 0 : chi = chi/sigma;
    indexPsiPi[j] = yMinMass + (j+0.5)*widthy;
    chiPsiPi[j] = chi;
    chi2PsiPi += chi*chi;
    if (chi != 0) dofPsiPi++;
  }

  dofPsiPi = dofPsiPi - (NORD1+1);
  probPsiPi = TMath::Prob(chi2PsiPi, dofPsiPi);
  
  TGraph* chiKPiGraph = new TGraph(xProj->GetNbinsX(), indexKPi, chiKPi);
  TGraph* chiPsiPiGraph = new TGraph(yProj->GetNbinsX(), indexPsiPi, chiPsiPi);

  //Angle efficiency fit

  int n1a = 0;
  if (erfSmooth == 1) n1a = 2;
  const int n2a = (NORD2+1)*(NORD3+1);
  double rectPar[2] = {6, 6}; 
  double parAng[n1a+n2a]; 

  if (erfSmooth == 1) 
    for (int index = 0; index < n1a; index++) 
      parAng[index] = rectPar[index];

  /*  if (readPar) {
      double parVal = 0;
      for(int i = n1a; i < n1a+n2a; i++){
      angParFile.good() ? angParFile >> parVal : parVal = 0;
      parAng[i] = parVal;
      }
      }
  */  else
    for (int index = n1a; index < n1a+n2a; index++)
      parAng[index] = 0;

  TF2* fitFunAng;
  if (erfSmooth == 1) fitFunAng = new TF2("fitAngEff", erfRectEff, xMinAng, xMaxAng, yMinAng, yMaxAng, n1a+n2a);
  else fitFunAng = new TF2("fitAngEff", rectEff, xMinAng, xMaxAng, yMinAng, yMaxAng, n2a);
  fitFunAng->SetParameters(parAng);

  for (int index = 0; index < n1a+n2a; index++)
    fitFun->SetParLimits(index, parAng[index]-margin, parAng[index]+margin);

  relEffTH2Ang->Fit("fitAngEff");

  //Angle efficiency Chi2D

  scalex = floor(nbins2/relEffTH2Ang->GetNbinsX());
  scaley = floor(nbins3/relEffTH2Ang->GetNbinsY());
  deltax = (xMaxAng - xMinAng)/nbins2;
  deltay = (yMaxAng - yMinAng)/nbins3;
  double widthxAng = (xMaxAng - xMinAng)/relEffTH2Ang->GetNbinsX();
  double widthyAng = (yMaxAng - yMinAng)/relEffTH2Ang->GetNbinsY();

  for (int i = 0; i < chiTH2Ang->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Ang->GetNbinsY(); j++) {
      sum = 0;
      for (int k = 0; k < scalex; k++) {
        for (int l = 0; l < scaley; l++) {
          xvar = xMinAng + i*widthxAng + (k + 0.5)*deltax;
          yvar = yMinAng + j*widthyAng + (l + 0.5)*deltay;
          sum += fitFunAng->Eval(xvar, yvar);
        }
      }
      sum *= pow(scalex*scaley, -1);
      chi = relEffTH2Ang->GetBinContent(i+1, j+1) - sum;
      sigma = relEffErrTH2Ang->GetBinContent(i+1, j+1);
      sigma == 0 ? chi = 0 : chi = pow(chi/sigma, 2);
      chiTH2Ang->SetBinContent(i+1, j+1, chi);

      if (chiThres && chi > thres) {
	chi = 0;
	relEffTH2AngCut->SetBinContent(i+1, j+1, 0);
	relEffErrTH2AngCut->SetBinContent(i+1, j+1, 0);
	nExcludedAng++;
      }
      chiTH2AngCut->SetBinContent(i+1, j+1, chi);
    }
  }

  chiAng = 0;
  dofAng = 0;
  for (int i = 0; i < chiTH2Ang->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Ang->GetNbinsY(); j++) {
      val = chiTH2AngCut->GetBinContent(i+1, j+1);
      chiAng += val;
      if (val != 0) dofAng++;
    }
  }

  dofAng = dofAng - (n1a+n2a);
  probAng = TMath::Prob(chiAng, dofAng);

  int nBinsHist1 = 48;
  TH1F* chi2AngDist = new TH1F("chi2AngDistribution", "chi2AngDistribution", nBinsHist1, 0, 12);
  for (int i = 0; i < chiTH2Ang->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2Ang->GetNbinsY(); j++) {
      chi2AngDist->Fill(chiTH2Ang->GetBinContent(i+1, j+1));
    }
  }
  chi2AngDist->AddBinContent(nBinsHist1, chi2AngDist->GetBinContent(nBinsHist1+1));
  TH1F* chi2AngDistCut = new TH1F("chi2AngDistributionCut", "chi2AngDistributionCut", nBinsHist1, 0, 12);
  for (int i = 0; i < chiTH2AngCut->GetNbinsX(); i++) {
    for (int j = 0; j < chiTH2AngCut->GetNbinsY(); j++) {
      chi2AngDistCut->Fill(chiTH2AngCut->GetBinContent(i+1, j+1));
    }
  }
  chi2AngDistCut->AddBinContent(nBinsHist1, chi2AngDistCut->GetBinContent(nBinsHist1+1));

  //Angles efficiency projections

  TH1F* xProjAng = new TH1F("cosMuMuEff", "cosMuMuEff", relEffTH2Ang->GetNbinsX(), xMinAng, xMaxAng);
  TH1F* yProjAng = new TH1F("phiEff", "phiEff", relEffTH2Ang->GetNbinsY(), yMinAng, yMaxAng);
  TH1F* xFitProjAng = new TH1F("cosMuMuEffFit","cosMuMuEffFit", nbins2, xMinAng, xMaxAng);
  TH1F* yFitProjAng = new TH1F("phiEffFit","phiEffFit", nbins3, yMinAng, yMaxAng);

  for (int i = 0; i < xProjAng->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffTH2Ang->GetNbinsY(); j++) {
      val = relEffTH2Ang->GetBinContent(i+1, j+1);
      sum += val;
    }
    xProjAng->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < yProjAng->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < relEffTH2Ang->GetNbinsX(); i++) {
      val = relEffTH2Ang->GetBinContent(i+1, j+1);
      sum += val;
    }
    yProjAng->SetBinContent(j+1, sum);
  }

  for (int i = 0; i < nbins2; i++) {
    sum = 0;
    xvar = xMinAng + (i + 0.5)*deltax;
    for (int j = 0; j < nbins3; j++) {
      yvar = yMinAng + (j + 0.5)*deltay;
      val = fitFunAng->Eval(xvar, yvar);
      sum += val;
    }
    sum *= pow(scaley, -1);
    xFitProjAng->SetBinContent(i+1, sum);
  }

  for (int j = 0; j < nbins3; j++) {
    sum = 0;
    yvar = yMinAng + (j + 0.5)*deltay;
    for (int i = 0; i < nbins2; i++) {
      xvar = xMinAng + (i + 0.5)*deltax;
      val = fitFunAng->Eval(xvar, yvar);
      sum += val;
    }
    sum *= pow(scalex, -1);
    yFitProjAng->SetBinContent(j+1, sum);
  }

  //Pull cosMuMu

  for (int i = 0; i < relEffErrTH2Ang->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < relEffErrTH2Ang->GetNbinsY(); j++) {
      sum += pow(relEffErrTH2Ang->GetBinContent(i+1, j+1), 2);
    }
    sum = pow(sum, 0.5);
    relEffErrCosMuMu->SetBinContent(i+1, sum);
  }

  chi2CosMuMu = 0;
  dofCosMuMu = 0;

  for (int i = 0; i < xProjAng->GetNbinsX(); i++) {
    sum = 0;
    for (int j = 0; j < yProjAng->GetNbinsX(); j++) {
      for (int k = 0; k < scalex; k++) {
        xvar = xMinAng + i*widthxAng + (k + 0.5)*deltax;
        for (int l = 0; l < scaley; l++) {
          yvar = yMinAng + j*widthyAng + (l + 0.5)*deltay;
          sum += fitFunAng->Eval(xvar, yvar);
        }
      }
    }
    sum *= pow(scalex*scaley, -1);
    chi = xProjAng->GetBinContent(i+1) - sum;
    sigma = relEffErrCosMuMu->GetBinContent(i+1);
    sigma == 0 ? chi = 0 : chi = chi/sigma;
    indexCosMuMu[i] = xMinAng + (i+0.5)*widthxAng;
    chiCosMuMu[i] = chi;
    chi2CosMuMu += chi*chi;
    if (chi != 0) dofCosMuMu++;
  }  

  dofCosMuMu = dofCosMuMu - (NORD2+1);
  probCosMuMu = TMath::Prob(chi2CosMuMu, dofCosMuMu);

  //Pull phi

  chi2Phi = 0;
  dofPhi = 0;

  for (int j = 0; j < relEffErrTH2Ang->GetNbinsY(); j++) {
    sum = 0;
    for (int i = 0; i < relEffErrTH2Ang->GetNbinsX(); i++) {
      sum += pow(relEffErrTH2Ang->GetBinContent(i+1, j+1), 2);
    }
    sum = pow(sum, 0.5);
    relEffErrPhi->SetBinContent(j+1, sum);
  }

  for (int j = 0; j < yProjAng->GetNbinsX(); j++) {
    sum = 0;
    for (int i = 0; i < xProjAng->GetNbinsX(); i++) {
      for (int k = 0; k < scalex; k++) {
        xvar = xMinAng + i*widthxAng + (k + 0.5)*deltax;
        for (int l = 0; l < scaley; l++) {
          yvar = yMinAng + j*widthyAng + (l + 0.5)*deltay;
          sum += fitFunAng->Eval(xvar, yvar);
        }
      }
    }
    sum *= pow(scalex*scaley, -1);
    chi = yProjAng->GetBinContent(j+1) - sum;
    sigma = relEffErrPhi->GetBinContent(j+1);
    sigma == 0 ? chi = 0 : chi = chi/sigma;
    indexPhi[j] = yMinAng + (j+0.5)*widthyAng;
    chiPhi[j] = chi;
    chi2Phi += chi*chi;
    if (chi != 0) dofPhi++;
  }
 
  dofPhi = dofPhi - (NORD3+1);
  probPhi = TMath::Prob(chi2Phi, dofPhi);

  TGraph* chiCosMuMuGraph = new TGraph(xProjAng->GetNbinsX(), indexCosMuMu, chiCosMuMu);
  TGraph* chiPhiGraph = new TGraph(yProjAng->GetNbinsX(), indexPhi, chiPhi);

  //Print chi-squared

  cout << "TH2Mass Chi2: " << chiMass << " - dof: " << dofMass << " - gof: " << probMass << endl;
  cout << "TH2Ang  Chi2: " << chiAng << " - dof: " << dofAng <<  " - gof: " << probAng << endl;
  //  cout << "TH1KPi  Chi2: " << chi2KPi << " - dof: " << dofKPi << " - gof: " << probKPi << endl;
  //  cout << "TH1PsiPi Chi2: " << chi2PsiPi << " - dof: " << dofPsiPi << " - gof: " << probPsiPi << endl;
  //  cout << "TH1CosMuMu Chi2: " << chi2CosMuMu << " - dof: " << dofCosMuMu << " - gof: " << probCosMuMu << endl;
  //  cout << "TH1Phi Chi2: " << chi2Phi << " - dof: " << dofPhi << " - gof: " << probPhi << endl;

  if (chiRoc) {
    outDalRoc << thres << "\t" << probMass << "\t" << nExcludedMass << endl;
    outAngRoc << thres << "\t" << probAng << "\t" << nExcludedAng << endl;
  }

  outDalRoc.close();
  outAngRoc.close();

  //Plot

  relEffTH2Mass->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  relEffTH2Mass->GetYaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  relEffTH2Ang->GetXaxis()->SetTitle("cos(MuMu)");
  relEffTH2Ang->GetYaxis()->SetTitle("phi [rad]");
  relEffTH2MassCut->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  relEffTH2MassCut->GetYaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  relEffTH2AngCut->GetXaxis()->SetTitle("cos(MuMu)");
  relEffTH2AngCut->GetYaxis()->SetTitle("phi [rad]");

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

  chiTH2Mass->SetTitle("#chi^{2} mass efficiency");
  chiTH2Mass->GetZaxis()->SetRangeUser(0, 10);
  chiTH2Ang->SetTitle("#chi^{2} angle efficiency");
  chiTH2Ang->GetZaxis()->SetRangeUser(0, 10);
  chiTH2Mass->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  chiTH2Mass->GetYaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  chiTH2Ang->GetXaxis()->SetTitle("cos(MuMu)");
  chiTH2Ang->GetYaxis()->SetTitle("phi [rad]");
  chiTH2MassCut->SetTitle("#chi^{2} mass efficiency");
  chiTH2MassCut->GetZaxis()->SetRangeUser(0, 10);
  chiTH2AngCut->SetTitle("#chi^{2} angle efficiency");
  chiTH2AngCut->GetZaxis()->SetRangeUser(0, 10);
  chiTH2MassCut->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  chiTH2MassCut->GetYaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  chiTH2AngCut->GetXaxis()->SetTitle("cos(MuMu)");
  chiTH2AngCut->GetYaxis()->SetTitle("phi [rad]");

  dalitzGraph->SetMarkerStyle(7);
  chiKPiGraph->SetMarkerStyle(7);
  chiPsiPiGraph->SetMarkerStyle(7);
  chiCosMuMuGraph->SetMarkerStyle(7);
  chiPhiGraph->SetMarkerStyle(7);
  chiKPiGraph->SetMinimum(-8);
  chiKPiGraph->SetMaximum(+8);
  chiPsiPiGraph->SetMinimum(-8);
  chiPsiPiGraph->SetMaximum(+8);
  chiCosMuMuGraph->SetMinimum(-8);
  chiCosMuMuGraph->SetMaximum(+8);
  chiPhiGraph->SetMinimum(-8);
  chiPhiGraph->SetMaximum(+8);
  chiKPiGraph->GetXaxis()->SetRangeUser(xProj->GetXaxis()->GetXmin(), xProj->GetXaxis()->GetXmax());
  chiPsiPiGraph->GetXaxis()->SetRangeUser(yProj->GetXaxis()->GetXmin(), yProj->GetXaxis()->GetXmax());
  chiCosMuMuGraph->GetXaxis()->SetRangeUser(xProjAng->GetXaxis()->GetXmin(), xProjAng->GetXaxis()->GetXmax()); 
  chiPhiGraph->GetXaxis()->SetRangeUser(yProjAng->GetXaxis()->GetXmin(), yProjAng->GetXaxis()->GetXmax());
  chiKPiGraph->SetTitle("Pull bin-by-bin");
  chiPsiPiGraph->SetTitle("Pull bin-by-bin");
  chiCosMuMuGraph->SetTitle("Pull bin-by-bin");
  chiPhiGraph->SetTitle("Pull bin-by-bin");
  chiKPiGraph->GetXaxis()->SetTitle("m(K#pi) [GeV/c2]");
  chiPsiPiGraph->GetXaxis()->SetTitle("m(J#psi#pi) [GeV/c2]");
  chiCosMuMuGraph->GetXaxis()->SetTitle("cos(MuMu)");
  chiPhiGraph->GetXaxis()->SetTitle("phi [rad]");

  chi2MassDist->GetXaxis()->SetTitle("#chi^{2}_{bin}");
  chi2MassDist->GetYaxis()->SetTitle("counts");
  chi2AngDist->GetXaxis()->SetTitle("#chi^{2}_{bin}");
  chi2AngDist->GetYaxis()->SetTitle("counts");
  chi2MassDistCut->SetFillColor(2);
  chi2AngDistCut->SetFillColor(2);

  TLine *line0up = new TLine(xMinMass, +3, xMaxMass, +3);
  TLine *line1up = new TLine(yMinMass, +3, yMaxMass, +3);
  TLine *line2up = new TLine(xMinAng, +3, xMaxAng, +3);
  TLine *line3up = new TLine(yMinAng, +3, yMaxAng, +3);
  TLine *line0do = new TLine(xMinMass, -3, xMaxMass, -3);
  TLine *line1do = new TLine(yMinMass, -3, yMaxMass, -3);
  TLine *line2do = new TLine(xMinAng, -3, xMaxAng, -3);
  TLine *line3do = new TLine(yMinAng, -3, yMaxAng, -3);

  line0up->SetLineStyle(9);
  line0up->SetLineColor(2);
  line1up->SetLineStyle(9);
  line1up->SetLineColor(2);
  line2up->SetLineStyle(9);
  line2up->SetLineColor(2);
  line3up->SetLineStyle(9);
  line3up->SetLineColor(2);
  line0do->SetLineStyle(9);
  line0do->SetLineColor(2);
  line1do->SetLineStyle(9);
  line1do->SetLineColor(2);
  line2do->SetLineStyle(9);
  line2do->SetLineColor(2);
  line3do->SetLineStyle(9);
  line3do->SetLineColor(2);

  TCanvas *c0 = new TCanvas("canvas", "canvas", 800, 600);
  c0->cd();
  relEffTH2MassCut->Draw("COLZ");
  dalitzGraph->Draw("Psame");

  TCanvas *c1 = new TCanvas("canvas1", "canvas1", 800, 600);
  c1->cd();
  relEffTH2AngCut->Draw("COLZ");

  TCanvas *cxMass = new TCanvas("canvasDpProjX", "canvasDpProjX", 800, 600);
  cxMass->Divide(1, 2);
  cxMass->cd(1);
  xProj->Draw();
  xFitProj->Draw("same");
  xProj->Draw("same");
  cxMass->cd(2);
  chiKPiGraph->Draw("AP");
  line0up->Draw();
  line0do->Draw();

  TCanvas *cyMass = new TCanvas("canvasDpProjY", "canvasDpProjY", 800, 600);  
  cyMass->Divide(1, 2);
  cyMass->cd(1);
  yProj->Draw();
  yFitProj->Draw("same");
  yProj->Draw("same");
  cyMass->cd(2);
  chiPsiPiGraph->Draw("AP");
  line1up->Draw();
  line1do->Draw();

  TCanvas *cxAng = new TCanvas("canvasAngProjX", "canvasAngProjX", 800, 600);
  cxAng->Divide(1, 2);
  cxAng->cd(1);
  xProjAng->Draw();
  xFitProjAng->Draw("same");
  xProjAng->Draw("same");
  cxAng->cd(2);
  chiCosMuMuGraph->Draw("AP");
  line2up->Draw();
  line2do->Draw();

  TCanvas *cyAng = new TCanvas("canvasAngProjY", "canvasAngProjY", 800, 600);
  cyAng->Divide(1, 2);
  cyAng->cd(1);
  yProjAng->Draw();
  yFitProjAng->Draw("same");
  yProjAng->Draw("same");
  cyAng->cd(2);
  chiPhiGraph->Draw("AP");
  line3up->Draw();
  line3do->Draw();

  TCanvas *chi2Mass = new TCanvas("chi2Mass", "chi2Mass", 800, 600);
  chi2Mass->cd();
  chiTH2MassCut->Draw("lego2");

  TCanvas *chi2Ang = new TCanvas("chi2Ang", "chi2Ang", 800, 600);
  chi2Ang->cd();
  chiTH2AngCut->Draw("lego2");

  TCanvas *chi2DistCanvas = new TCanvas("chi2distribution", "chi2distribution", 800, 600);
  chi2DistCanvas->Divide(1, 2);
  chi2DistCanvas->cd(1);
  chi2DistCanvas->cd(1)->SetLogy();
  chi2MassDist->Draw();
  chi2MassDistCut->Draw("same");
  chi2DistCanvas->cd(2);
  chi2DistCanvas->cd(2)->SetLogy();
  chi2AngDist->Draw();
 chi2AngDistCut->Draw("same");

  return 0;
}
