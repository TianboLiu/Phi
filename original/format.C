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
    double W, q, Q, cth, t, t0, ds, stat, syst;
    ifstream file1("clasdb_E63M2.txt");//charged mode
    for (int i = 0; i < 8; i++) file1.getline(tmp, 256);
    FILE * f1000 = fopen("datasets/1000.csv","w");
    fprintf(f1000, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file1 >> W >> cth >> ds >> stat){
      if (W == 2.735 || W == 2.745) continue;
      n++;
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      syst = ds * 0.1115;
      fprintf(f1000, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dcth", ds, stat, syst, -syst,"mub");
    }
    fclose(f1000);
    cout << "1000.csv written" << endl;
    file1.close();
    ifstream file2("clasdb_E63M3.txt");//neutral mode
    for (int i = 0; i < 8; i++) file2.getline(tmp, 256);
    FILE * f1001 = fopen("datasets/1001.csv","w");
    fprintf(f1001, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    n = 0;
    while (file2 >> W >> cth >> ds >> stat){
      n++;
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q * cth;
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      syst = ds * 0.1135;
      fprintf(f1001, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dcth", ds, stat, syst, -syst,"mub");
    }
    fclose(f1001);
    cout << "1001.csv written" << endl;
    file2.close();
  }

  if (opt == 2){//zeus 1996
    double W, q, Q, cth, t, t0, ds, stat, syst, temp;
    ifstream file("zeus1996.txt");
    for (int i = 0; i < 8; i++) file.getline(tmp, 256);
    FILE * f1002 = fopen("datasets/1002.csv", "w");
    fprintf(f1002, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    W = 70.0;
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    int n = 0;
    while (file >> t >> temp >> temp >> ds >> stat >> temp){
      n++;
      t = -t;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      syst = sqrt(0.06*0.06+0.09*0.09) * ds;
      fprintf(f1002, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst,"mub/GeV2");
    }
    fclose(f1002);
    cout << "1002.csv written" << endl;
    file.close();
  }

  if (opt == 3){//zeus 2000
    double W, q, Q, cth, t, t0, ds, stat, syst1, syst2, syst3, syst4, temp;
    ifstream file("zeus2000.txt");
    for (int i = 0; i < 9; i++) file.getline(tmp, 256);
    FILE * f1003 = fopen("../datasets/1003.csv", "w");
    fprintf(f1003, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    W = 94.0;
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    int n = 0;
    while (file >> t >> temp >> temp >> ds >> stat >> temp >> syst1 >> syst2 >> syst3 >> syst4){
      n++;
      t = -t;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1003, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, 
	      sqrt(syst1*syst1+syst3*syst3+0.15*0.15*ds*ds), -sqrt(syst2*syst2+syst4*syst4+0.15*0.15*ds*ds),
	      "mub/GeV2");
    }
    fclose(f1003);
    cout << "1003.csv written" << endl;
    file.close();
  }

  if (opt == 4){//desy094
    double Eg[8] = {3.0, 3.5, 4.0, 4.5, 5.1, 5.6, 6.2, 6.7};
    double W[7], q, Q, t, t0, cth, ds[7], stat[7], temp;
    ifstream file("desy094.txt");
    for (int i = 0; i < 8; i++) file.getline(tmp, 256);
    FILE * f1004 = fopen("../datasets/1004.csv", "w");
    fprintf(f1004, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    for (int i = 0; i < 7; i++) 
      W[i] = 0.5 * (sqrt(pow(Eg[i] + Mp, 2) - Eg[i] * Eg[i]) + sqrt(pow(Eg[i+1] + Mp, 2) - Eg[i+1] * Eg[i+1]));
    int n = 0;
    while (file >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp >> ds[2] >> stat[2] >> temp >> ds[3] >> stat[3] >> temp >> ds[4] >> stat[4] >> temp >> ds[5] >> stat[5] >> temp >> ds[6] >> stat[6] >> temp){
      for (int i = 0; i < 7; i++){
	if (ds[i] > 0){
	  n++;
	  t = -t;
	  q = (W[i] * W[i] - Mp * Mp) / (2.0 * W[i]);
	  Q = sqrt((W[i] * W[i] - pow(Mp + Mphi, 2)) * (W[i] * W[i] - pow(Mp - Mphi, 2))) / (2.0 * W[i]);
	  t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
	  cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
	  fprintf(f1004, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
		  n, W[i], q, Q, cth, t, t0, "ds/dt", ds[i], stat[i], 0.04 * ds[i], -0.04 * ds[i], "mub/GeV2"); 
	}
      }
    }
    fclose(f1004);
    cout << "1004.csv written" << endl;
    file.close();
  }

  if (opt == 5){//CERN-WA-004
    ifstream file("cern-wa-004.txt");
    for (int i = 0; i < 7; i++) file.getline(tmp, 256);
    double Eg, q, Q, stat, syst, sig, temp;
    file >> Eg >> temp >> temp >> sig >> stat >> temp >> syst >> temp;
    FILE * f1005 = fopen("../datasets/1005.csv", "w");
    fprintf(f1005, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    double W = sqrt(pow(Eg + Mp, 2) - Eg * Eg);
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    double t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    fprintf(f1005, "%d,%.6E,%.6E,%.6E,%s,%s,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	    1, W, q, Q, "NA", "NA", t0, "sig*Br", sig/1000, stat/1000, syst/1000, -syst/1000, "mub");
    fclose(f1005);
    cout << "1005.csv written" << endl;
    file.close();
  }

  if (opt == 6){//BONN 1974
    double W, q, Q, cth, t, t0, ds, stat, syst, temp;
    ifstream file("bonn1974a.txt");
    for (int i = 0; i < 9; i++) file.getline(tmp, 256);
    FILE * f1006 = fopen("../datasets/1006.csv", "w");
    fprintf(f1006, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    W = sqrt(pow(2.0+Mp, 2) - 2.0 * 2.0);
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    int n = 0;
    while (file >> temp >> t >> temp >> stat >> syst >> ds >> temp){
      n++;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      ds = ds / 1000.0;
      fprintf(f1006, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds, ds*stat/100, ds*syst/100, -ds*syst/100, "mub/GeV2");
    }
    fclose(f1006);
    cout << "1006.csv written" << endl;
    file.close();
  }

  if (opt == 7){//LAMP2
    double W, q, Q, cth, t, t0, ds[3], stat[3], temp;
    ifstream file1("lamp2a.txt");
    for (int i = 0; i < 9; i++) file1.getline(tmp, 256);
    FILE * f1007 = fopen("../datasets/1007.csv", "w");
    fprintf(f1007, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file1 >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp){
      n++;
      t = -t;
      W = 0.5 * (2.476 + 2.83);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[0], stat[0], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (3.144 + 2.83);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[1], stat[1], 0., -0., "mub/GeV2");
    }
    file1.close();
    ifstream file2("lamp2b.txt");
    for (int i = 0; i < 9; i++) file2.getline(tmp, 256);
    while (file1 >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp >> ds[2] >> stat[2] >> temp){
      n++;
      t = -t;
      W = 0.5 * (2.476 + 2.695);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[0], stat[0], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (2.695 + 2.896);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[1], stat[1], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (2.896 + 3.144);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[2], stat[2], 0., -0., "mub/GeV2");
    }
    file2.close();
    fclose(f1007);
    cout << "1007.csv written" << endl;
  }
      
	     

  return 0;
}
