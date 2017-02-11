#include "mcmc.h"

int main(){

  char filename[] = "datafit.dat";
  CreateEmptySet(filename);
  //AddSet(1000, filename);//CLAS
  //AddSet(1001, filename);//CLAS
  AddSet(1002, filename);
  AddSet(1003, filename);
  AddSet(1004, filename);
  AddSet(1006, filename);
  AddSet(1007, filename);
  AddSet(1008, filename);
  AddSet(1010, filename);
  AddSet(1012, filename);
  //AddSet(1013, filename);//LEPS
  //AddSet(1015, filename);//CLAS
  AddSet(1016, filename);
  AddSet(1017, filename);
  AddSet(1018, filename);
  AddSet(1019, filename);
  AddSet(1020, filename);
  AddSet(1021, filename);

  LoadData(filename, 1);

  double par[20] = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01, 0., 0., 0., 0.,
		    1.01, 1.01, 0.01, 0.01, 0.01, 0., 0., 0., 0., 0.};
  //double par[6] = {3.5650E-01,3.5566E+00,7.2062E-01,8.1056E-03,2.4265E+00,2.4562E-03};
  Model = &Model1;

  cout << Minimizer(par, 20) << endl;
  //MarkovChainScan(par, 6);
  //MetropolicScan(par, 6, 1000);

  ClearAll();
  return 0;
}
