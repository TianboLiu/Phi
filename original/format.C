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
    FILE * f1003 = fopen("datasets/1003.csv", "w");
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
    FILE * f1004 = fopen("datasets/1004.csv", "w");
    fprintf(f1004, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    for (int i = 0; i < 7; i++) 
      W[i] = 0.5 * (sqrt(pow(Eg[i] + Mp, 2) - Eg[i] * Eg[i]) + sqrt(pow(Eg[i+1] + Mp, 2) - Eg[i+1] * Eg[i+1]));
    int n = 0;
    while (file >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp >> ds[2] >> stat[2] >> temp >> ds[3] >> stat[3] >> temp >> ds[4] >> stat[4] >> temp >> ds[5] >> stat[5] >> temp >> ds[6] >> stat[6] >> temp){
      t = -t;
      for (int i = 0; i < 7; i++){
	if (ds[i] > 0){
	  n++;
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
    FILE * f1005 = fopen("datasets/1005.csv", "w");
    fprintf(f1005, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    double W = 6.197;
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
    FILE * f1006 = fopen("datasets/1006.csv", "w");
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
    FILE * f1007 = fopen("datasets/1007.csv", "w");
    fprintf(f1007, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file1 >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp){
      n++;
      t = -t;
      W = 0.5 * (2.476 + 2.83);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      if (t > -0.05) t = 0.5 * (-0.05 + t0);
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[0], stat[0], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (3.144 + 2.83);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      if (t > -0.05) t = 0.5 * (-0.05 + t0);
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[1], stat[1], 0., -0., "mub/GeV2");
    }
    file1.close();
    ifstream file2("lamp2b.txt");
    for (int i = 0; i < 9; i++) file2.getline(tmp, 256);
    while (file2 >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp >> ds[2] >> stat[2] >> temp){
      n++;
      t = -t;
      W = 0.5 * (2.476 + 2.695);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      if (t > -0.05) t = 0.5 * (-0.05 + t0);
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[0], stat[0], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (2.695 + 2.896);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      if (t > -0.05) t = 0.5 * (-0.05 + t0);
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[1], stat[1], 0., -0., "mub/GeV2");
      n++;
      W = 0.5 * (2.896 + 3.144);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      if (t > -0.05) t = 0.5 * (-0.05 + t0);
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1007, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      n, W, q, Q, cth, t, t0, "ds/dt", ds[2], stat[2], 0., -0., "mub/GeV2");
    }
    file2.close();
    fclose(f1007);
    cout << "1007.csv written" << endl;
  }

  if (opt == 8){//SLAC 1973
    double W, q, Q, cth, t, t0, ds[2], stat[2], temp;
    ifstream file("slac1973.txt");
    for (int i = 0; i < 8; i++) file.getline(tmp, 256);
    FILE * f1008 = fopen("datasets/1008.csv", "w");
    fprintf(f1008, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file >> t >> temp >> temp >> ds[0] >> stat[0] >> temp >> ds[1] >> stat[1] >> temp){
      t = -t;
      W = sqrt(pow((2.8+4.7)/2.0 + Mp, 2) - pow((2.8+4.7)/2.0, 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      if (ds[0] > 0)
	fprintf(f1008, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds[0], stat[0], 0.0, -0.0, "mub/GeV2");
      W = sqrt(pow(9.3 + Mp, 2) - pow(9.3, 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      if (ds[1] > 0)
	fprintf(f1008, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds[1], stat[1], 0.0, -0.0, "mub/GeV2");
    }
    fclose(f1008);
    cout << "1008.csv written" << endl;
    file.close();
  }

  if (opt == 9){//Fermilab E401
    double Eg, W, q, Q, cth, t, t0, ds, stat, temp;
    ifstream file1("fermilab-e401a.txt");
    for (int i = 0; i < 3; i++) file1.getline(tmp, 256);
    FILE * f1009 = fopen("datasets/1009.csv", "w");
    fprintf(f1009, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file1 >> Eg >> temp >> temp >> ds >> stat){
      W = sqrt(pow(Eg + Mp, 2) - Eg * Eg);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      fprintf(f1009, "%d,%.6E,%.6E,%.6E,%s,%s,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, "NA", "NA", t0, "sig", ds, stat, 0.08*ds, -0.08*ds, "mub");
    }
    fclose(f1009);
    cout << "1009.csv written" <<  endl;
    file1.close();
    ifstream file2("fermilab-e401b.txt");
    for (int i = 0; i < 3; i++) file2.getline(tmp, 256);
    FILE * f1010 = fopen("datasets/1010.csv", "w");
    fprintf(f1009, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    n = 0;
    W = sqrt(pow(100.0 + Mp, 2) - 100.0 * 100.0);
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    while (file2 >> t >> temp >> temp >> ds >> stat){
      t = -t;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      fprintf(f1010, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, 0.08*ds, -0.08*ds, "mub/GeV2");
    }
    fclose(f1010);
    cout << "1010.csv written" << endl;
    file2.close();
  }

  if (opt == 10){//Fermilab E25
    double Eg, W, q, Q, t0, ds, stat, temp;
    ifstream file("fermilab-e25.txt");
    for (int i = 0; i < 7; i++) file.getline(tmp, 256);
    FILE * f1011 = fopen("datasets/1011.csv", "w");
    fprintf(f1011, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    while (file >> Eg >> temp >> temp >> temp >> temp >> temp >> ds >> stat >> temp){
      W = sqrt(pow(Eg + Mp, 2) - Eg * Eg);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      fprintf(f1011, "%d,%.6E,%.6E,%.6E,%s,%s,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, "NA", "NA", t0, "sig", ds, stat, 0.05*ds, -0.05*ds, "mub");
    }
    fclose(f1011);
    cout << "1011.csv written" <<  endl;
    file.close();
  }

  if (opt == 11){//CLAS E93-031
    double W, q, Q, cth, t, t0, ds, stat;
    ifstream file("clasdb_E77M1.txt");
    for (int i = 0; i < 8; i++) file.getline(tmp, 256);
    FILE * f1012 = fopen("datasets/1012.csv", "w");
    fprintf(f1012, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    W = sqrt(pow(3.6 + Mp, 2) - 3.6 * 3.6);
    while (file >> t >> ds >> stat){
      t = -t;
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      ds = ds / 1000.0;
      stat = stat / 1000.0;
      fprintf(f1012, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, 0.0, -0.0, "mub/GeV2");
    }
    fclose(f1012);
    cout << "1012.csv written" << endl;
    file.close();
  }

  if (opt == 12){//LEPS
    double W, q, Q, cth, t, t0, ds, stat, syst, temp, W1, W2;
    string files[8] = {"leps2005a.txt", "leps2005b.txt", "leps2005c.txt", "leps2005d.txt", "leps2005e.txt", "leps2005f.txt", "leps2005g.txt", "leps2005h.txt"};
    FILE * f1013 = fopen("datasets/1013.csv", "w");
    fprintf(f1013, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    syst = sqrt(0.8 * 0.8 + 2.1 * 2.1 + 3.0 * 3.0) / 100.0;
    int n = 0;
    for (int j = 0; j < 8; j++){
      ifstream file(files[j].data());
      for (int i = 0; i < 5; i++) file.getline(tmp, 256);
      file >> tmp >> tmp >> W1 >> tmp >> W2 >> tmp >> tmp;
      W = 0.5 * (W1 + W2);
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      for (int i = 0; i < 3; i++) file.getline(tmp, 256);
      while (file >> t >> temp >> temp >> ds >> stat >> temp){
	t = t + t0;
	cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
	fprintf(f1013, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
		++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, ds*syst, -ds*syst, "mub/GeV2");
      }
      file.close();
    }
    fclose(f1013);
    cout << "1013.csv written" << endl;
  }
     
  if (opt == 13){//CERN-WA-057
    ifstream file("cern-wa-057.txt");
    for (int i = 0; i < 7; i++) file.getline(tmp, 256);
    double Eg, q, Q, stat, syst, sig, temp;
    file >> Eg >> temp >> temp >> sig >> stat >> temp >> syst >> temp;
    FILE * f1014 = fopen("datasets/1014.csv", "w");
    fprintf(f1014, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    double W = 6.197;
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    double t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    fprintf(f1014, "%d,%.6E,%.6E,%.6E,%s,%s,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	    1, W, q, Q, "NA", "NA", t0, "sig*Br", sig/1000, stat/1000, syst/1000, -syst/1000, "mub");
    fclose(f1014);
    cout << "1014.csv written" << endl;
    file.close();
  }

  if (opt == 14){//CLAS neutral
    string files[9] = {"clas-e04-021a.txt", "clas-e04-021b.txt", "clas-e04-021c.txt", "clas-e04-021d.txt", "clas-e04-021e.txt", "clas-e04-021f.txt", "clas-e04-021g.txt", "clas-e04-021h.txt", "clas-e04-021i.txt"};
    double Eg, W, q, Q, cth, tt[5], t[5], t0, ds[5], syst[5], stat[5];
    FILE * f1015 = fopen("datasets/1015.csv", "w");
    fprintf(f1015, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    for (int j = 0; j < 9; j++){
      ifstream file(files[j].data());
      file >> tt[0] >> tt[1] >> tt[2] >> tt[3] >> tt[4];
      while (file >> Eg >> W >> ds[0] >> syst[0] >> stat[0] >> ds[1] >> syst[1] >> stat[1] >> ds[2] >> syst[2] >> stat[2] >> ds[3] >> syst[3] >> stat[3] >> ds[4] >> syst[4] >> stat[4]){
	q = (W * W - Mp * Mp) / (2.0 * W);
	Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
	t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
	for (int i = 0; i < 5; i++) t[i] = t0 - tt[i];
	for (int i = 0; i < 5; i++){
	  if (ds[i] > 0){
	    cth = (t[i] + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
	    fprintf(f1015, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
		    ++n, W, q, Q, cth, t[i], t0, "ds/dt", ds[i], stat[i], syst[i], -syst[i], "mub/GeV2");
	  }
	}
      }
      file.close();
    }
    fclose(f1015);
    cout << "1015.csv written" << endl;
  }

  if (opt == 15){//SAPHIR
    double E1, E2, W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    string files[4] = {"saphira.txt", "saphirb.txt", "saphirc.txt", "saphird.txt"};
    FILE * f1016 = fopen("datasets/1016.csv", "w");
    fprintf(f1016, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    for (int j = 0; j < 4; j++){
      ifstream file(files[j].data());
      file >> tmp >> E1 >> E2 >> tmp;
      file.getline(tmp, 256);
      file.getline(tmp, 256);
      W = sqrt(pow(0.5 * (E1 + E2) + Mp, 2) - pow(0.5 * (E1 + E2), 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      while (file >> t >> ds){
	file >> temp >> dsup;
	file >> temp >> dsdown;
	t = t0 - t;
	cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
	stat = 0.5 * (dsup - dsdown);
	syst = 0.1 * ds;
	fprintf(f1016, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
		++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
      }
      file.close();
    }
    fclose(f1016);
    cout << "1016.csv written" << endl;
  }

  if (opt == 16){//Cornell 1972
    double W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    ifstream file("cornell1972.txt");
    for (int i = 0; i < 3; i++) file.getline(tmp, 256);
    FILE * f1017 = fopen("datasets/1017.csv", "w");
    fprintf(f1017, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    W = sqrt(pow(8.5 + Mp, 2) - 8.5 * 8.5);
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    while (file >> t >> ds){
      file >> temp >> dsup;
      file >> temp >> dsdown;
      t = -t;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      stat = 0.5 * (dsup - dsdown);
      syst = 0.0;
      fprintf(f1017, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
              ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
    }
    fclose(f1017);
    cout << "1017.csv written" << endl;
    file.close();
  }

  if (opt == 17){//ABBHHM (DESY HBC)
    double E1, E2, W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    string files[2] = {"abbhhma.txt", "abbhhmb.txt"};
    FILE * f1018 = fopen("datasets/1018.csv", "w");
    fprintf(f1018, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    for (int j = 0; j < 2; j++){
      ifstream file(files[j].data());
      file >> tmp >> E1 >> E2 >> tmp;
      file.getline(tmp, 256);
      file.getline(tmp, 256);
      W = sqrt(pow(0.5 * (E1 + E2) + Mp, 2) - pow(0.5 * (E1 + E2), 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      while (file >> t >> ds){
        file >> temp >> dsup;
        file >> temp >> dsdown;
        t = - t;
        cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
        stat = 0.5 * (dsup - dsdown);
        syst = 0.0;
        fprintf(f1018, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
                ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
      }
      file.close();
    }
    fclose(f1018);
    cout << "1018.csv written" << endl;
  }

  if (opt == 18){//Cornell 1971
    double E1, W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    string files[3] = {"cornell1971a.txt", "cornell1971b.txt", "cornell1971c.txt"};
    FILE * f1019 = fopen("datasets/1019.csv", "w");
    fprintf(f1019, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    for (int j = 0; j < 3; j++){
      ifstream file(files[j].data());
      file >> tmp >> E1 >> tmp;
      file.getline(tmp, 256);
      file.getline(tmp, 256);
      W = sqrt(pow(E1 + Mp, 2) - pow(E1, 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      while (file >> t >> ds){
        file >> temp >> dsup;
        file >> temp >> dsdown;
        t = - t;
        cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
        stat = 0.5 * (dsup - dsdown);
        syst = 0.1 * ds;
        fprintf(f1019, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
                ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
      }
      file.close();
    }
    fclose(f1019);
    cout << "1019.csv written" << endl;
  }

  if (opt == 19){//SLAC 1970
    double E1, W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    string files[7] = {"slac1970a.txt", "slac1970b.txt", "slac1970c.txt", "slac1970d.txt", "slac1970e.txt", "slac1970f.txt", "slac1970g.txt"};
    FILE * f1020 = fopen("datasets/1020.csv", "w");
    fprintf(f1020, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    for (int j = 0; j < 7; j++){
      ifstream file(files[j].data());
      file >> tmp >> E1 >> tmp;
      file.getline(tmp, 256);
      file.getline(tmp, 256);
      W = sqrt(pow(E1 + Mp, 2) - pow(E1, 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      while (file >> t >> ds){
        file >> temp >> dsup;
        file >> temp >> dsdown;
        t = - t;
        cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
        stat = 0.5 * (dsup - dsdown);
        syst = 0.15 * ds;
        fprintf(f1020, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
                ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
      }
      file.close();
    }
    fclose(f1020);
    cout << "1020.csv written" << endl;
  }

  if (opt == 20){//SLAC 1973L
    double E1, W, q, Q, cth, t, t0, ds, dsup, dsdown, stat, syst, temp;
    FILE * f1021 = fopen("datasets/1021.csv", "w");
    fprintf(f1021, "i,W,q,Q,cth,t,t0,obs,value,stat,syst+,syst-,unit\n");
    int n = 0;
    ifstream file1("slac1973La.txt");
    file1 >> tmp >> t >> tmp;
    file1.getline(tmp, 256);
    file1.getline(tmp, 256);
    t = -t;
    while (file1 >> E1 >> ds){
      file1 >> temp >> dsup;
      file1 >> temp >> dsdown;
      W = sqrt(pow(E1 + Mp, 2) - pow(E1, 2));
      q = (W * W - Mp * Mp) / (2.0 * W);
      Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
      t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      stat = 0.5 * (dsup - dsdown);
      syst = 0.0;
      fprintf(f1021, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
    }
    file1.close();
    ifstream file2("slac1973Lb.txt");
    file2 >> tmp >> E1 >> tmp;
    file2.getline(tmp, 256);
    file2.getline(tmp, 256);
    W = sqrt(pow(E1 + Mp, 2) - pow(E1, 2));
    q = (W * W - Mp * Mp) / (2.0 * W);
    Q = sqrt((W * W - pow(Mp + Mphi, 2)) * (W * W - pow(Mp - Mphi, 2))) / (2.0 * W);
    t0 = Mphi * Mphi - 2.0 * q * sqrt(Mphi * Mphi + Q * Q) + 2.0 * q * Q;
    while (file2 >> t >> ds){
      file2 >> temp >> dsup;
      file2 >> temp >> dsdown;
      t = -t;
      cth = (t + 2.0 * q * sqrt(Mphi * Mphi + Q * Q) - Mphi * Mphi) / (2.0 * q * Q);
      stat = 0.5 * (dsup - dsdown);
      syst = 0.0;
      fprintf(f1021, "%d,%.6E,%.6E,%.6E,%.6E,%.6E,%.6E,%s,%.6E,%.6E,%.6E,%.6E,%s\n",
	      ++n, W, q, Q, cth, t, t0, "ds/dt", ds, stat, syst, -syst, "mub/GeV2");
    }
    file2.close();
    fclose(f1021);
    cout << "1021.csv written" << endl;
  }

  return 0;
}
