#include "hmc.h"

int main(const int argc, const char * argv[]){

  DataSet * dataset = new DataSet();
  dataset->AddSet(1001);
  int obs;
  double W, t, ds, err;
  dataset->GetData(12, obs, W, t, ds, err, 2);
  cout << obs << " " << W << " " << t << " "<< ds << " " << err << endl;

  HMC * hmc = new HMC();
  hmc->SetModel(&Model0);

  return 0;
}
