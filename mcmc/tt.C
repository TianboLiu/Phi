#include "mcmc.h"

int main(){

  // CreateEmptySet("datafit.dat");
  // AddSet(1002, "datafit.dat");
  // AddSet(1003, "datafit.dat");

  LoadData("datafit.dat", 0);

  double par[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};;
  Model = &Model0;

  cout << Minimizer(par, 2) << endl;
  
  ClearAll();
  return 0;
}
