#ifndef _MCMC_H_
#define _MCMC_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <omp.h>

#include "TROOT.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TMatrixDEigen.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

const double Mp = 0.938272;//proton mass
const double Mphi = 1.019455;//phi meson mass
const double W0 = Mp + Mphi;//minimal invariant mass

int readline(ifstream * file, string * varstr, const int col);
int CreateEmptySet(const char * filename);
int AddSet(const int ID, const char * filename);
int LoadData(const char * filename, const int errtype);
double Model0(const double * var, const double * par);
double Model1(const double * var, const double * par);
double Model2(const double * var, const double * par);
double (*Model)(const double *, const double *);
double Chi2(const double * par);
double Minimizer(const double * par0, const int Npar, double * res, const char * minName, const char * algoName);
double MetropolicScan(const double * par0, const int Npar, const int steps);
double MarkovChainScan(const double * par0, const int Npar);
int CheckParameters(const double * par, const int Npar);
double IndividualSetChi2(const int ID, const int errtype, const double * par);
void ClearAll();
double * _W, * _t, * _cth, * _ds, * _err, * _var2, * _model, * _chi2;
int * _ID, * _obs;
int _Npoints;



int readline(ifstream * file, string * varstr, const int col = 13){//read line from csv data file
  char tmp[300];
  for (int i = 0; i < col - 1; i++){
    file->getline(tmp, 50, ',');
    varstr[i] = tmp;
  }
  if (file->good()){
    file->getline(tmp, 50, '\n');
    varstr[col-1] = tmp;
    return 1;
  }
  else
    return 0;
}

int CreateEmptySet(const char * filename){//create data set with no entries
  FILE * fp = fopen(filename, "w");
  fclose(fp);
  return 1;
}

int AddSet(const int ID, const char * filename){//add data set to file for fitting
  ostringstream IDstr;
  IDstr << ID;
  string fdata = "datasets/" + IDstr.str() + ".csv";
  ifstream fread(fdata.data());
  if (!fread.is_open()){
    printf("Data set %d is not found in datasets!\n", ID);
    return 0;
  }
  char tmp[300];
  FILE * fadd = fopen(filename, "a");
  fread.getline(tmp, 300);
  string varstr[13];
  while (readline(&fread, varstr, 13)){
    fprintf(fadd, "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
	    ID, 
	    varstr[1].data(), varstr[2].data(), varstr[3].data(), varstr[4].data(), varstr[5].data(), varstr[6].data(), varstr[7].data(), varstr[8].data(), varstr[9].data(), varstr[10].data(), varstr[11].data(), varstr[12].data());
  }
  fread.close();
  fclose(fadd);
  return 1;
}

int LoadData(const char * filename, const int errtype){//load data
  ifstream file(filename);
  if (!file.is_open()){
    printf("Fail to file datafile for fit %s!\n", filename);
    return 0;
  }
  _Npoints = 0;
  char tmp[300];
  while (file.getline(tmp, 300)){
    _Npoints++;
  }
  int ID, obs;
  double W, cth, t, ds, stat, syst1, syst2, err, temp;
  char name[10];
  _ID = new int [_Npoints];
  _obs = new int [_Npoints];
  _W = new double [_Npoints];
  _cth = new double [_Npoints];
  _t = new double [_Npoints];
  _ds = new double [_Npoints];
  _err = new double [_Npoints];
  _var2 = new double [2 * _Npoints];
  _model = new double [_Npoints];
  _chi2 = new double [_Npoints];
  file.clear();
  file.seekg(0, ios::beg);
  printf("Loading data from %s ...\n", filename);
  for (int i = 0; i < _Npoints; i++){
    file >> ID >> W >> temp >> temp >> cth >> t >> temp >> name >> ds >> stat >> syst1 >> syst2 >> tmp;
    if (strcmp(name, "ds/dt") == 0) obs = 0;
    else if (strcmp(name, "ds/dcth") == 0) obs = 1;
    else obs = -1;
    if (errtype == 0) err = stat;
    else if (errtype == 1) err = sqrt(pow(stat, 2) + 0.5 * pow(syst1, 2) + 0.5 * pow(syst2, 2));
    else err = -1.;
    _ID[i] = ID;
    _obs[i] = obs;
    _W[i] = W;
    _cth[i] = cth;
    _t[i] = t;
    _ds[i] = ds;
    _err[i] = err;
    _var2[2*i] = W;
    _var2[2*i+1] = t;
  }
  printf("%d data loaded\n", _Npoints);
  file.close();
  return 1;
}

double Chi2(const double * par){
#pragma omp parallel num_threads(8)
  {
    #pragma omp for schedule(dynamic, 3)
    for (int i = 0; i < _Npoints; i++){
      _model[i] = Model(&_var2[2*i], par);
      if (_obs[i] == 1){
	_model[i] = _model[i] * 2.0 * (_W[i] * _W[i] - Mp * Mp) / (2.0 * _W[i])
	  * sqrt((_W[i] * _W[i] - pow(Mp + Mphi, 2)) * (_W[i] * _W[i] - pow(Mp - Mphi, 2))) / (2.0 * _W[i]);
      }
      _chi2[i] = pow( (_model[i] - _ds[i]) / _err[i], 2);
    }
  }
  double sum = 0.0;
  for (int i = 0; i < _Npoints; i++)
    sum += _chi2[i];
  return sum;
}

double Minimizer(const double * par0, const int Npar, double * res, const char * minName = "Minuit2", const char * algoName = "Migrad"){//least chi-square search
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(10000);
  min->SetTolerance(1.0e-3);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2, Npar);
  min->SetFunction(f);
  string parName;
  ostringstream pari;
  for (int i = 0; i < Npar; i++){
    pari.str("");
    pari << i;
    parName = "a" + pari.str();
    min->SetVariable(i, parName.data(), par0[i], 1.0e-5);
    if (par0[i] == 0.0 || par0[i] == 1.0)
      min->FixVariable(i);
  }
  min->Minimize();
  const double * X = min->X();
  for (int i = 0; i < Npar; i++){
    res[i] = X[i];
  }
  return min->MinValue() / (_Npoints - Npar);
}

int CheckParameters(const double * par, const int Npar){
  for (int i = 0; i < Npar; i++)
    if (par[i] < 0) return 0;
  return 1;
}

double Model0(const double * var, const double * par){
  double W = var[0];
  double t = var[1];
  double q = (W * W - Mp * Mp) / (2.0 * W);
  double Q = sqrt((W * W - pow(Mp + Mphi, 2)) *  (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
  double t0 = -2.0 * q * sqrt(Mphi * Mphi + Q * Q) + Mphi * Mphi + 2.0 * q * Q;
  double Wr = (W - W0) / W0;
  double tr = (t0 - t) / (W0 * W0);
  double Aterm = atan(par[1] * pow(Wr, par[2])) * (1.0 + par[3] * Wr / (pow(Wr - par[4], 2) + par[5] * par[5]) + par[6] * log(1.0 + Wr));
  double Bterm = exp(-par[7] * tr);
  double ds = par[0] * Aterm * Bterm;
  return ds;
}

double Model1(const double * var, const double * par){
  double a[10], b[10];
  for (int i = 0; i < 10; i++){
    a[i] = par[i];
    b[i] = par[i+10];
  }
  double W = var[0];
  double t = var[1];
  double q = (W * W - Mp * Mp) / (2.0 * W);
  double Q = sqrt((W * W - pow(Mp + Mphi, 2)) *  (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
  double t0 = -2.0 * q * sqrt(Mphi * Mphi + Q * Q) + Mphi * Mphi + 2.0 * q * Q;
  double sr = (W * W - W0 * W0) / (Mp * Mp);
  double tr = (t0 - t) / (Mp * Mp);
  double A0 = a[0] * atan(a[1] * pow(sr, a[2])) * (1.0 + a[3] / (pow(sr - a[4],2) + pow(a[4] * a[5],2)));
  double B0 = b[0] * pow(sr, b[1]);
  double B1 = 1.0 + b[3] * tr + b[4] * tr * tr;
  double ds = 100 * A0 * exp(-B0 * pow(tr, 1.0 - b[2])) * B1;
  return ds;
}

double Model2(const double * var, const double * par){
  return 2.0;
}

double IndividualSetChi2(const int ID, const int errtype, const double * par){
  CreateEmptySet("oneset.dat");
  AddSet(ID, "oneset.dat");
  LoadData("oneset.dat", errtype);
  double sum = Chi2(par);
  return sum / _Npoints;
}


void ClearAll(){//release memory
  delete[] _ID;
  delete[] _obs;
  delete[] _W;
  delete[] _cth;
  delete[] _t;
  delete[] _ds;
  delete[] _err;
  delete[] _var2;
  delete[] _model;
  delete[] _chi2;
}

#endif
