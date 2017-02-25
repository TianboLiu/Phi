#ifndef _HMC_H_
#define _HMC_H_

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
#include "Math/Derivator.h"
#include "TMatrixDEigen.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

using namespace std;

const double Mp = 0.938272;//proton mass
const double Mphi = 1.019455;//phi meson mass
const double W0 = Mp + Mphi;//minimal invariant mass

class DataSet{
 private:
  int Npoint;
  int * index;
  int * obs;
  double * var[3];
  double * value;
  double * err0;
  double * err1;
  double * err2;
  TString dir;
  TString stmp;
 public:
  DataSet();
  int AddSet(const int setid);
  int PrintData(const int n) const;
  int GetData(const int n, int& obs, double& W, double& t, double& ds, double& err, const int etype) const;
  int GetNpoint() const;
};

DataSet::DataSet(){
  Npoint = 0;
  const int Nmax = 4000;
  dir = "datasets/";
  index = new int[Nmax];
  obs = new int[Nmax];
  var[0] = new double[Nmax];
  var[1] = new double[Nmax];
  var[2] = new double[Nmax];
  value = new double[Nmax];
  err0 = new double[Nmax];
  err1 = new double[Nmax];
  err2 = new double[Nmax];
}

int DataSet::AddSet(const int setid){//add a data set
  TString stmp = dir + Form("%d", setid) + ".csv";//data file name
  ifstream infile(stmp.Data());//open data file
  if (!infile.is_open()){
    printf("File %s does not exist!\n", stmp.Data());
    return -1;
  }
  stmp.Clear();
  stmp.ReadFile(infile);//read file into stmp
  infile.close();//close data file
  stmp.ReplaceAll(",", " ");//replace , with space
  FILE * ftmp = fopen("stmp.dat", "w");//open a tmp dat file
  stmp.Puts(ftmp);//write data into tmp dat file
  fclose(ftmp);//save tmp dat file
  ifstream file("stmp.dat");//open tmp dat file
  TString sobs, st;
  double W, ds, stat, syst1, syst2;
  stmp.ReadLine(file);//skip the header line in file
  while (file >> stmp >> W >> stmp >> stmp >> stmp >> st >> stmp >> sobs >> ds >> stat >> syst1 >> syst2 >> stmp){
    if (!sobs.EqualTo("ds/dt") && !sobs.EqualTo("ds/dcth")){
      printf("Set %d is not a differential cross section data!\n", setid);
      file.close();
      return 0;
    }
    index[Npoint] = setid;
    if (sobs.EqualTo("ds/dt"))
      obs[Npoint] = 0;
    else
      obs[Npoint] = 1;
    var[0][Npoint] = W;
    var[1][Npoint] = st.Atof();
    value[Npoint] = ds;
    err0[Npoint] = stat;
    err1[Npoint] = 0.5 * (syst1 - syst2);
    err2[Npoint] = sqrt(pow(err0[Npoint], 2) + pow(err1[Npoint], 2));
    Npoint++;
  }
  printf("Set %d added.\n", setid);
  file.close();
  return 0;
}
  
int DataSet::GetData(const int n, int& opt, double& W, double& t, double& ds, double& err, const int etype = 2) const{
  if (!(n < Npoint)){
    perror("Overflow dataset size!");
    return -1;
  }
  opt = obs[n];
  W = var[0][n];
  t = var[1][n];
  ds = value[n];
  if (etype == 0)
    err = err0[n];
  else if (etype == 1)
    err = err1[n];
  else if (etype == 2)
    err = err2[n];
  else {
    perror("Invalid etype!");
    return -1;
  }
  return 0;
}

int DataSet::GetNpoint() const{
  return Npoint;
}

//////////////////////////////////////////////////////

class HMC{
 private:
  DataSet * dataset;
  int dim;
  double step;
  int length;
  double * Q;
  double * P;
 public:
  HMC();
  int SetDataSet(DataSet * data);
  double (*Model)(const double *, const double *);
  int SetModel(double (*f)(const double *, const double *));
  double Hamiltonian(const double * q, const double * p);
  double Potential(const double * q);
  int GradPotential(const double * q, double * grad);
  double Kinetic(const double * p);
};
HMC hh;

HMC::HMC(){
  dim = 20;
  Q = new double[dim];
  P = new double[dim];
}

int HMC::SetDataSet(DataSet * data){
  dataset = data;
  return 0;
}

int HMC::SetModel(double (*f)(const double * var, const double *)){
  Model = f;
  return 0;
}

double HMC::Potential(const double * q){
  return 0;
}

int HMC::GradPotential(const double * q, double * grad){
  ROOT::Math::Functor f(this, &HMC::Potential, dim);
  //ROOT::Math::Derivator df;
  //df.SetFunction(f);
  for (int i = 0; i < dim; i++){
    // grad[i] = df.Eval(&f, q, i);
  }
  return 0;
}



//////////////////////////////////////////////////////

double Model0(const double * var, const double * par){
  return 0;
}

  
#endif
