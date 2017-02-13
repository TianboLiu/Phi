#include "mcmc.h"

int main(int argc, char * argv[]){

  char filename[] = "datafit.dat";
  CreateEmptySet(filename);
  AddSet(1000, filename);//CLAS
  AddSet(1001, filename);//CLAS
  AddSet(1002, filename);
  AddSet(1003, filename);
  AddSet(1004, filename);
  AddSet(1006, filename);
  AddSet(1007, filename);
  AddSet(1008, filename);
  AddSet(1010, filename);
  AddSet(1012, filename);
  AddSet(1013, filename);//LEPS
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

  int fix1[20] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
		 0, 0, 0, 0, 0, 0, 0, 1, 1, 1};
  double par1[20] = {0.043, 0.190, 0.44, 0.001, 3.07, 0.0314, 0., 0., 0., 0.,
		     2.004, 0.243, -0.0313, -0.942, 0.644, 0.1, 0.5, 0., 0., 0.};

  double res[20];
  Model = &Model1;

  cout << Minimizer(par1, fix1, 20, res) << endl;
  //MarkovChainScan(par, 6);
  //MetropolicScan(par, 6, 1000);
  //res[6] = 0.01;
  int fix2[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		 1, 1, 1, 1, 1, 0, 0, 1, 1, 1};;
  //cout << Minimizer(res, fix2, 20, par) << endl;

  if (argc > 1){
    int ni = atoi(argv[1]);
    cout << ni << " " << IndividualSetChi2(ni, 1, res) << endl;
  }

  ClearAll();
  return 0;
}
