#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

int main(int argc, char * argv[]){
  
  if (argc < 2){
    printf("./script <opt> <datafile>\n");
    printf("opt = 1: clasdb_E63");
    return 1;
  }

  int opt = atoi(argv[1]);
  char tmp[256];
  
  const double Mp = 0.938272;
  const double Mphi = 1.019455;

  if (opt == 1){
    double W, q, Q, cth, t, ds, stat, syst;
    ifstream file1("clasdb_E63M2.txt");//charged mode
    for (int i = 0; i < 8; i++) file1.getline(tmp, 256);
    FILE * f1000 = fopen("../datasets/1000.csv","w");
    fprintf(f1000, "W,q,Q,cth,t,ds,stat,syst,obs,unit,experiment\n");
    while (file1 >> W >> cth >> ds >> stat){
      if (W == 2.735 || W == 2.745) continue;
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      syst = ds * 0.1115;
      fprintf(f1000, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, syst,
	      "ds/dcth", "mb", "clas2013");
    }
    fclose(f1000);
    file1.close();
    ifstream file2("clasdb_E63M3.txt");//neutral mode
    for (int i = 0; i < 8; i++) file2.getline(tmp, 256);
    FILE * f1001 = fopen("../datasets/1001.csv","w");
    fprintf(f1001, "W,q,Q,cth,t,ds,stat,syst,obs,unit,experiment\n");
    while (file1 >> W >> cth >> ds >> stat){
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      syst = ds * 0.1135;
      fprintf(f1001, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, syst,
	      "ds/dcth", "mb", "clas2013");
    }
    fclose(f1001);
    file2.close();
  }

  return 0;
}
