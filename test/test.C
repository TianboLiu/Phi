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

  return 0;
}
