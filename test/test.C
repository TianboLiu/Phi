#include "test.h"


int main(int argc, char * argv[]){
  if (argc < 2){
    cout << "./test <opt>" << endl;
    return 1;
  }

  int opt = atoi(argv[1]);

  if (opt == 1){//
    MakeData(100);
  }

  if (opt == 2){//
    printdata();
  }

  if (opt == 3){//
    plotdata();
  }

  if (opt == 4){//
    mcmc1();
  }

  if (opt == 5){//
    plotchain();
  }

  if (opt == 6){//
    minimizer();
  }


    

  return 0;
}
