#include "mcmc.h"

int main(){

  char filename[] = "datafit.dat";
  CreateEmptySet(filename);
  AddSet(1000, filename);
  //AddSet(1001, filename);
  //AddSet(1015, filename);
  //AddSet(1002, filename);
  //AddSet(1003, filename);

  LoadData(filename, 1);

  double par[10] = {1.0, 1.01, 1.01, 1.01, 0.5, 0.01, 0.01, 0.01, 0.01, 0.01};;
  //double par[6] = {3.5650E-01,3.5566E+00,7.2062E-01,8.1056E-03,2.4265E+00,2.4562E-03};
  Model = &Model1;

  cout << Minimizer(par, 7) << endl;
  //MarkovChainScan(par, 6);
  //MetropolicScan(par, 6, 1000);

  ClearAll();
  return 0;
}
