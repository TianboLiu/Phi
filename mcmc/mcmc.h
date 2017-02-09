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
double Minimizer(const double * par0, const int Npar, const char * minName, const char * algoName);
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

double Minimizer(const double * par0, const int Npar, const char * minName = "Minuit2", const char * algoName = "Migrad"){//least chi-square search
  ROOT::Math::Minimizer * min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(1.0e-6);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&Chi2, Npar);
  min->SetFunction(f);
  string parName;
  ostringstream pari;
  for (int i = 0; i < Npar; i++){
    pari.str("");
    pari << i;
    parName = "a" + pari.str();
    min->SetVariable(i, parName.data(), par0[i], 1.0e-4);
  }
  min->Minimize();
  return min->MinValue() / (_Npoints - Npar);
}

double Model0(const double * var, const double * par){
  double W = var[0];
  double t = var[1];
  double q = (W * W - Mp * Mp) / (2.0 * W);
  double Q = sqrt((W * W - pow(Mp + Mphi, 2)) *  (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
  double t0 = -2.0 * q * sqrt(Mphi * Mphi + Q * Q) + Mphi * Mphi + 2.0 * q * Q;
  double ds = par[0] * exp(par[1] * (t - t0));
  return ds;
}

double Model1(const double * var, const double * par){
  return 1.0;
}

double Model2(const double * var, const double * par){
  return 2.0;
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
