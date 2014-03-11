#include <iostream>
#include <cmath>
#include <complex>
#include <new>
#include <memory>
#include <algorithm>

// #include <fstream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"

// #include <boost/algorithm/string.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>
// #include <boost/lexical_cast.hpp>

using namespace std;

typedef std::complex<double> cplx;

int main(int argc, char const *argv[])
{
  matrixReal M(5,6);
  M.makeIdentity();
  M(2,1) = 10.0;

  matrixReal N(6,2);
  N.makeIdentity();
  N(4,1) = -3.0;
  N(2,0) = -20.0;

  cout << M << N;

  matrixReal O = M*N;
  cout << O << O.transpose();

  matrixReal Q = N-N*2.0;
  cout << Q;

  matrixReal R(8,8);
  R.makeIdentity();
  R(3,4) = -2.0;
  R(4,3) = -2.0;
  R(0,4) = 3.0;
  R(4,0) = 3.0;
  cout << R;
  cout << R.getSub(1,1,4,4);


  vector <double> vals(8,0.0);
  R.diagonalize(vals.data());
  cout << R;
  vector <double> svd_vals(8,0.0);
  tuple<shared_ptr<matrixReal>,shared_ptr<matrixReal>> S = R.svd(svd_vals);
  matrixReal T = *get<0>(S);
  matrixReal U = *get<1>(S);
  cout << T << U;
  T.setSub(1,0,O.transpose());
  cout << T;
  return 0;
}
