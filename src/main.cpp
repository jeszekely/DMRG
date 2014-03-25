#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>

#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <memory>


#include "matrix.hpp"
#include "block.hpp"
#include "dmrg.hpp"
#include "davidson.hpp"
// #include <boost/algorithm/string.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>
// #include <boost/lexical_cast.hpp>

using namespace std;

typedef std::complex<double> cplx;

int main(int argc, char const *argv[])
{
  int maxChainLen=20; //This means the entire thing is 20, so 10 on each side
  vector<int> maxKeepNum = {10, 20, 30};

  infiniteSystem infChain(maxChainLen, maxKeepNum);

  //infChain.runInfinite();
  
  finiteSystem finiteChain(maxChainLen, maxKeepNum);
  finiteChain.sweep();

  return 0;
}
