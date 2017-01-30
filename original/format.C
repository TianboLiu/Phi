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

  if (opt == 1){//clas 2013
    double W, q, Q, cth, t, ds, stat, syst;
    ifstream file1("clasdb_E63M2.txt");//charged mode
    for (int i = 0; i < 8; i++) file1.getline(tmp, 256);
    FILE * f1000 = fopen("../datasets/1000.csv","w");
    fprintf(f1000, "W,q,Q,cth,t,ds,stat,syst+,syst-,obs,unit,experiment\n");
    while (file1 >> W >> cth >> ds >> stat){
      if (W == 2.735 || W == 2.745) continue;
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      syst = ds * 0.1115;
      fprintf(f1000, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, syst, syst,
	      "ds/dcth", "mb", "clas2013");
    }
    fclose(f1000);
    file1.close();
    ifstream file2("clasdb_E63M3.txt");//neutral mode
    for (int i = 0; i < 8; i++) file2.getline(tmp, 256);
    FILE * f1001 = fopen("../datasets/1001.csv","w");
    fprintf(f1001, "W,q,Q,cth,t,ds,stat,syst+,syst-,obs,unit,experiment\n");
    while (file2 >> W >> cth >> ds >> stat){
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      syst = ds * 0.1135;
      fprintf(f1001, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, syst, syst,
	      "ds/dcth", "mb", "clas2013");
    }
    fclose(f1001);
    file2.close();
  }

  if (opt == 2){//zeus 1996
    double W, q, Q, cth, t, ds, stat, syst;
    ifstream file("zeus1996.txt");
    for (int i = 0; i < 3; i++) file.getline(tmp, 256);
    FILE * f1002 = fopen("../datasets/1002.csv", "w");
    fprintf(f1002, "W,q,Q,cth,t,ds,stat,syst+,syst-,obs,unit,experiment\n");
    W = 70.0;
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    while (file >> t >> ds >> stat){
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      syst = sqrt(0.06*0.06+0.09*0.09) * ds;
      fprintf(f1002, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, syst, syst,
	      "ds/dt", "mb/GeV2", "zeus1996");
    }
    fclose(f1002);
    file.close();
  }

  if (opt == 3){//zeus 2000
    double W, q, Q, cth, t, ds, stat, syst1, syst2, syst3, syst4;
    ifstream file("zeus2000.txt");
    for (int i = 0; i < 3; i++) file.getline(tmp, 256);
    FILE * f1003 = fopen("../datasets/1003.csv", "w");
    fprintf(f1003, "W,q,Q,cth,t,ds,stat,syst+,syst-,obs,unit,experiment\n");
    W = 95.0;
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    while (file >> t >> ds >> stat >> syst1 >> syst2 >> syst3 >> syst4){
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1003, "%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%s,%s\n",
	      W, q, Q, cth, t, ds, stat, 
	      sqrt(syst1*syst1 + syst3*syst3 + 0.15*0.15*ds*ds), 
	      sqrt(syst2*syst2 + syst4*syst4 + 0.15*0.15*ds*ds),
	      "ds/dt", "mb/GeV2", "zeus2000");
    }
    fclose(f1003);
    file.close();
  }


  return 0;
}
