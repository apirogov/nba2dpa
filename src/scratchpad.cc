#include <iostream>
#include "io.hh"
using namespace std;
using namespace nbautils;

//this is for trying stuff out. do whatever you want here.
int main(int argc, char *argv[]) {
  string file="";
  if (argc>1)
    file = argv[1];
  auto bas = parse_hoa_ba(file);
  if (bas.empty())
    cout << "something went wrong!" << endl;

  auto &ba = bas.front();
  if (argc>2){
    auto ps = powerset_product(*ba);
    print_hoa(*ps);
  } else {
    auto ps = powerset_construction(*ba);
    print_hoa(*ps);
  }

}
