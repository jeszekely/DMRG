#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "block.hpp"
#include "dmrg.hpp"
#include "davidson.hpp"
// #include <boost/algorithm/string.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>
// #include <boost/lexical_cast.hpp>

using namespace std;

int main(int argc, char const *argv[])
{

#if 0
  matrixReal M(5,6);
  M.makeIdentity();
  M(2,1) = 10.0;
  M.printMem();
  matrixReal N(6,2);
  N.makeIdentity();
  N(4,1) = -3.0;
  N(2,0) = -20.0;

  cout << M << N;

  matrixReal O = M*N;
  cout << O << *O.transpose();

  matrixReal Q = N-N*2.0;
  cout << Q;

  matrixReal R(8,8);
  R.makeIdentity();
  R(3,4) = -2.0;
  R(4,3) = -2.0;
  R(0,4) = 3.0;
  R(4,0) = 3.0;
  cout << R;
  cout << *R.getSub(1,1,4,4);
  R.printMem();
  vector <double> vals(8,0.0);
  R.diagonalize(vals.data());
  cout << "Diagonalized Matrix: " << endl << R << endl;
  cout << "Eigenvalues :" << endl << *vals.data() << endl;

  vector <double> svd_vals(8,0.0);
  shared_ptr<matrixReal> T, U;
  tie(T,U) = R.svd(svd_vals);
  cout << *T << *U;
  T->setSub(1,0,*O.transpose());
  cout << *T;

  vector <double> HDiags(U->nc(),0.0);
  for (int ii = 0; ii < U->nc(); ii++)
  HDiags[ii] = U->element(ii,ii);
  Davidson Diag(2,8,4,HDiags,[&](matrixReal& A){return R*A;});
  vector<double> EVals(8);
  Diag.diagonalize(EVals);
#endif

#if 1
  int maxChainLen=20; //This means the entire thing is 20, so 10 on each side
  vector<int> maxKeepNum = {10, 20, 30};
  finiteSystem finiteChain(maxChainLen, maxKeepNum);
  finiteChain.sweep();
#endif

  return 0;
}
