#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "TRandom3.h"
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

int mcmc1(){
  loaddata();
  TRandom3 * random = new TRandom3(0);
  double par0[2] = {2.0, 2.0};
  double par1[2];
  FILE * fc = fopen("chain1.txt", "w");
  double xx0, P0, xx1, P1;
  xx0 = chi2(par0);
  P0 = exp(-0.5 * xx0);
  for (int n = 0; n < 1000; n++){
    fprintf(fc, "%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    printf("%d\t %.6E\t %.6E %.6E %.6E\n",
	    n, par0[0], par0[1], xx0, P0);
    while (true){
      par1[0] = par0[0] + random->Gaus(0.0, 10.0);
      par1[1] = par0[1] + random->Gaus(0.0, 10.0);
      xx1 = chi2(par1);
      P1 = exp(-0.5 * xx1);
      if (random->Uniform(0.0, 1.0) < P1 / P0) break;
    }
    par0[0] = par1[0];
    par0[1] = par1[1];
    xx0 = xx1;
    P0 = P1;
  }
  fclose(fc);
  return 0;
}

int plotchain(){
  ifstream fchain("chain1.txt");
  double t, a, b, xx, P;
  int n = 0;
  TGraph * achain = new TGraph(1000);
  TGraph * bchain = new TGraph(1000);
  while (fchain >> t >> a >> b >> xx >> P){
    achain->SetPoint(n, t, a);
    bchain->SetPoint(n, t, b);
    n++;
  }
  fchain.close();
  achain->SetLineColor(4);
  bchain->SetLineColor(2);
  TCanvas * c1 = new TCanvas("c1", "", 800, 600);
  achain->Draw("aL");
  bchain->Draw("Lsame");
  c1->Print("c1.pdf");
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
  TGraphErrors * gd = new TGraphErrors(100, _x, _y, 0, _err);
  gd->SetMarkerStyle(8);
  gd->SetMarkerSize(0.4);
  gd->SetLineWidth(2);
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
    err = sqrt(pow(f0(&x, par0up) - f0(&x, par0down), 2) + pow(f0(&x, par1up) - f0(&x, par1down), 2)) / 2.0;
    gfcov->SetPointError(i, 0, err);
  }
  gfcov->SetFillColor(4);
  TCanvas * c2 = new TCanvas("c2", "", 800, 600);
  gd->Draw("ape");
  gfcov->Draw("4same");
  c2->Print("c2.pdf");
  return 0;
}
		    
