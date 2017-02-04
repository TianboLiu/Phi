#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "TRandom3.h"

int MakeData(const int NN = 100){
  TRandom3 * random = new TRandom3(1);

  FILE * file = fopen("pseudodata.csv", "w");
  double x, y, err;
  for (int i = 0; i < NN; i++){
    x = random->Uniform(1.0, 10.0);
    err = random->Uniform(0.1, 0.5);
    y = random->Gaus(x, err);
    fprintf(file, "%d,%.6E,%.6E,%.6E\n",
	    i, x, y, err);
  }
  fclose(file);
  cout << "pseudodata.csv written, with " << NN << " points" << endl;
  return 0;
}
    
