#ifndef _TEST_H_
#define _TEST_H_

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "TRandom3.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMatrixDEigen.h"
#include "TH1D.h"
#include "TF1.h"

int MakeData(const int NN = 100){
  TRandom3 * random = new TRandom3(1);

  FILE * file = fopen("pseudodata.csv", "w");
  double x, y, err;
  for (int i = 0; i < NN; i++){
    x = random->Uniform(1.0, 5.0);
    y = exp(-x);
    err = random->Uniform(0.1, 0.5) * y;
    y = y + random->Gaus(0, err);
    fprintf(file, "%d,%.6E,%.6E,%.6E\n",
	    i, x, y, err);
  }
  fclose(file);
  cout << "pseudodata.csv written, with " << NN << " points" << endl;
  return 0;
}

int readline(ifstream * fp, double * var){
  char tmp[100];
  fp->getline(tmp, 100, ',');
  var[0] = atof(tmp);
  fp->getline(tmp, 100, ',');
  var[1] = atof(tmp);
  fp->getline(tmp, 100, ',');
  var[2] = atof(tmp);
  if (!fp->good()) return 0;
  fp->getline(tmp, 100, '\n');
  var[3] = atof(tmp);
  return 1;
}

int printdata(){
  double var[4];
  ifstream file("pseudodata.csv");
  while (readline(&file, var)){
    cout << var[0] << " " << var[1] << " " << var[2] << " " << var[3] << endl;
  }
  file.close();
  return 0;
}

int plotdata(){
  double x[100], y[100], err[100];
  double var[4];
  ifstream file("pseudodata.csv");
  for (int i = 0; i < 100; i++){
    if (readline(&file, var)){
      x[i] = var[1];
      y[i] = var[2];
      err[i] = var[3];
    }
  }
  TGraphErrors * gd = new TGraphErrors(100, x, y, 0, err);
  gd->SetMarkerStyle(8);
  gd->SetLineWidth(1.5);
  TCanvas * c0 = new TCanvas("c0", "", 800, 600);
  gd->Draw("ape");
  c0->Print("c0.pdf");
  return 0;
}

double f0(const double * x, const double * par){
  return par[0] * exp(-par[1] * x[0]);
}

double f1(const double * x, const double * par){
  return par[0] * exp(-(par[1] + par[2]) * x[0]);
}

double _x[100], _y[100], _err[100]; 
int loaddata(){
  double var[4];
  ifstream file("pseudodata.csv");
  for (int i = 0; i < 100; i++){
    if (readline(&file, var)){
      _x[i] = var[1];
      _y[i] = var[2];
      _err[i] = var[3];
    }
  }
  file.close();
  return 0;
}

double chi2(const double * par){
  double sum = 0.0;
  for (int i = 0; i < 100; i++){
    sum += pow((f0(&_x[i], par) - _y[i]) / _err[i], 2);
  }
  return sum;
}

int mcmc0(const double da, const double db, const char * savename = "chain0.txt"){
  loaddata();
  TRandom3 * random = new TRandom3(0);
  double par0[2] = {2.0, 0.0};
  double par1[2];
  FILE * fc = fopen(savename, "w");
  double xx0, P0, xx1, P1;
  xx0 = chi2(par0);
  P0 = exp(-0.5 * xx0);
  for (int n = 0; n < 10000; n++){
    fprintf(fc, "%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    printf("%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    par1[0] = par0[0] + random->Uniform(-da, da);
    par1[1] = par0[1] + random->Uniform(-db, db);
    xx1 = chi2(par1);
    P1 = exp(-0.5 * xx1);
    if (random->Uniform(0.0, 1.0) < exp(-0.5 * (xx1 - xx0))){
      par0[0] = par1[0];
      par0[1] = par1[1];
      xx0 = xx1;
      P0 = P1;
    }
  }
  fclose(fc);
  return 0;
}

int mcmc1(const double sa, const double sb, const char * savename = "chain1.txt"){
  loaddata();
  TRandom3 * random = new TRandom3(0);
  double par0[2] = {2.0, 0.0};
  double par1[2];
  FILE * fc = fopen(savename, "w");
  double xx0, P0, xx1, P1;
  xx0 = chi2(par0);
  P0 = exp(-0.5 * xx0);
  for (int n = 0; n < 10000; n++){
    fprintf(fc, "%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    printf("%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    par1[0] = par0[0] + random->Gaus(0.0, sa);
    par1[1] = par0[1] + random->Gaus(0.0, sb);
    xx1 = chi2(par1);
    P1 = exp(-0.5 * xx1);
    if (random->Uniform(0.0, 1.0) < exp(-0.5 * (xx1 - xx0))){
      par0[0] = par1[0];
      par0[1] = par1[1];
      xx0 = xx1;
      P0 = P1;
    }
  }
  fclose(fc);
  return 0;
}

int plotchain(const char * chainfile, const char * savefile){
  ifstream fchain(chainfile);
  double t, a, b, xx, P;
  int n = 0;
  TGraph * achain = new TGraph(10000);
  TGraph * bchain = new TGraph(10000);
  TH1D * afull = new TH1D("afull", "", 60, 0.7, 1.3);
  afull->SetDirectory(0);
  TH1D * apost = new TH1D("apost", "", 60, 0.7, 1.3);
  TH1D * bfull = new TH1D("bfull", "", 60, 0.85, 1.15);
  TH1D * bpost = new TH1D("bpost", "", 60, 0.85, 1.15);
  while (fchain >> t >> a >> b >> xx >> P){
    achain->SetPoint(n, t, a);
    bchain->SetPoint(n, t, b);
    n++;
    afull->Fill(a);
    bfull->Fill(b);
    if (t >= 2000){
      apost->Fill(a);
      bpost->Fill(b);
    }
  }
  fchain.close();
  achain->SetLineColor(4);
  bchain->SetLineColor(2);
  afull->SetLineColor(4);
  apost->SetLineColor(4);
  apost->SetLineWidth(2);
  bfull->SetLineColor(2);
  bpost->SetLineColor(2);
  bpost->SetLineWidth(2);
  gStyle->SetOptStat(0);
  TCanvas * c1 = new TCanvas("c1", "", 1600, 1200);
  c1->Divide(2,2);
  c1->cd(1);
  achain->GetXaxis()->SetRangeUser(0,2000);
  achain->DrawClone("aL");
  bchain->DrawClone("Lsame");
  c1->cd(2);
  achain->GetXaxis()->SetRangeUser(0,10000);
  achain->DrawClone("aL");
  bchain->DrawClone("Lsame");
  c1->cd(3);
  afull->DrawClone();
  apost->DrawClone("same");
  c1->cd(4);
  bfull->DrawClone();
  bpost->DrawClone("same");
  c1->Print(savefile);
  afull->Delete();
  apost->Delete();
  bfull->Delete();
  bpost->Delete();
  achain->Delete();
  bchain->Delete();
  c1->Close();
  return 0;
}

int comparisonplotnaive(const char * chainfile){
  loaddata();
  TH1D * base = new TH1D("base", "", 1, 0.0, 5.0);
  base->GetYaxis()->SetRangeUser(-0.1, 1.1);
  TGraphErrors * gd = new TGraphErrors(100, _x, _y, 0, _err);
  gd->SetMarkerStyle(8);
  gd->SetMarkerSize(0.8);
  gd->SetLineWidth(1);
  gd->SetLineColor(1);
  TGraphErrors * gmc = new TGraphErrors(500);
  double x, y, error, par[2];
  double t, a, b, xx, P;
  ifstream fc(chainfile);
  TH1D * hc = new TH1D("hc", "", 3, 0.0, 3.0);
  for (int i = 0; i < 500; i++){
    x = 0.01 * i;
    fc.clear();
    fc.seekg(0, ios::beg);
    while (fc >> t >> a >> b >> xx >> P){
      if (t >= 2000){
	par[0] = a;
	par[1] = b;
	hc->Fill(0.5, f0(&x, par));
      }
    }
    y = hc->GetBinContent(1) / 8000.0;
    gmc->SetPoint(i, x, y);
    fc.clear();
    fc.seekg(0, ios::beg);
    while (fc >> t >> a >> b >> xx >> P){
      if (t >= 2000){
	par[0] = a;
	par[1] = b;
	hc->Fill(1.5, pow(f0(&x, par) - y, 2));
      }
    }
    error = sqrt(hc->GetBinContent(2) / (8000.0 - 1.0));
    gmc->SetPointError(i, 0, error);
  }
  fc.close();
  gmc->SetFillColor(2);
  TCanvas * c3 = new TCanvas("c3", "", 800, 600);
  base->Draw();
  gd->Draw("pesame");
  gmc->Draw("4same");
  c3->Print("c3.pdf");
  return 0;
}

int minimizer(const char * minName = "Minuit2", const char * algoName = "Migrad"){
  loaddata();
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&chi2, 2);
  min->SetFunction(f);
  min->SetVariable(0, "a", 1.0, 1.0e-4);
  min->SetVariable(1, "b", 1.0, 1.0e-4);
  min->Minimize();
  const double * xs = min->X();
  double par[2] = {xs[0], xs[1]};
  cout << par[0] << "\t" << par[1] << "\t" << chi2(par) << "\t" << exp(-0.5 * chi2(par)) << endl;
  double cov[4];
  min->GetCovMatrix(cov);
  TMatrixD Cov(2,2);
  Cov(0,0) = cov[0]; Cov(0,1) = cov[1]; Cov(1,0) = cov[2]; Cov(1,1) = cov[3];
  TMatrixDEigen Eigen(Cov);
  Cov.Print();
  TMatrixD Values = Eigen.GetEigenValues();
  TMatrixD Vecs = Eigen.GetEigenVectors();
  Values.Print();
  Vecs.Print();
  return 0;
}

int comparisonplot(){
  loaddata();
  TH1D * base = new TH1D("base", "", 1, 0.0, 5.0);
  base->GetYaxis()->SetRangeUser(-0.1, 1.1);
  TGraphErrors * gd = new TGraphErrors(100, _x, _y, 0, _err);
  gd->SetMarkerStyle(8);
  gd->SetMarkerSize(0.8);
  gd->SetLineWidth(1);
  gd->SetLineColor(1);
  char tmp[100];
  ifstream fmini("cov.txt");
  fmini.getline(tmp, 100);
  double minipar[2] = {1.09403, 1.03312};
  double l0 = sqrt(0.004493);
  double l1 = sqrt(3.712e-5);
  double dpar0[2] = {0.9683, 0.2496};
  double dpar1[2] = {-0.2496, 0.9683};
  double par0up[2] = {minipar[0] + l0 * dpar0[0], minipar[1] + l0 * dpar0[1]};
  double par0down[2] = {minipar[0] - l0 * dpar0[0], minipar[1] - l0 * dpar0[1]};
  double par1up[2] = {minipar[0] + l1 * dpar1[0], minipar[1] + l1 * dpar1[1]};
  double par1down[2] = {minipar[0] - l1 * dpar1[0], minipar[1] - l1 * dpar1[1]};
  TGraphErrors * gfcov = new TGraphErrors(500);
  double x, error;
  for (int i = 0; i < 500; i++){
    x = 0.01 * i;
    gfcov->SetPoint(i, x, f0(&x, minipar));
    error = sqrt(pow(f0(&x, par0up) - f0(&x, par0down), 2) + pow(f0(&x, par1up) - f0(&x, par1down), 2)) / 2.0;
    gfcov->SetPointError(i, 0, error);
  }
  gfcov->SetFillColor(4);
  TCanvas * c2 = new TCanvas("c2", "", 800, 600);
  base->Draw();
  gd->Draw("pesame");
  gfcov->Draw("4same");
  c2->Print("c2.pdf");
  return 0;
}


#endif
