#include "test.h"


int main(int argc, char * argv[]){
  if (argc < 2){
    cout << "./test <opt>" << endl;
    return 1;
  }

  int opt = atoi(argv[1]);

  if (opt == 1){//
    //MakeData(100);
  }

  if (opt == 2){//
    printdata();
  }

  if (opt == 3){//
    plotdata();
  }

  if (opt == 4){//
    mcmc0(1.0, 1.0, "chain0_1.txt");
    mcmc0(0.1, 0.1, "chain0_0.1.txt");
    mcmc0(0.01, 0.01, "chain0_0.01.txt");
    mcmc1(1.0, 1.0, "chain1_1.txt");
    mcmc1(0.1, 0.1, "chain1_0.1.txt");
    mcmc1(0.01, 0.01, "chain1_0.01.txt");
  }

  if (opt == 5){//
    plotchain("chain0_1.txt","chain0_1.pdf");
    plotchain("chain0_0.1.txt","chain0_0.1.pdf");
    plotchain("chain0_0.01.txt","chain0_0.01.pdf");
    plotchain("chain1_1.txt","chain1_1.pdf");
    plotchain("chain1_0.1.txt","chain1_0.1.pdf");
    plotchain("chain1_0.01.txt","chain1_0.01.pdf");
  }

  if (opt == 6){//
    minimizer();
  }

  if (opt == 7){//
    comparisonplot();
  }

  if (opt == 8){//
    comparisonplotnaive("chain0_0.1.txt");
  }


    

  return 0;
}
