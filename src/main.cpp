#include <iostream>
#include <cmath>
#include <complex>
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
#include "dmrg.hpp"

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
  cout << R.getSub(1,1,4,4);
  R.printMem();
  vector <double> vals(8,0.0);
  R.diagonalize(vals.data());
  cout << R;
  vector <double> svd_vals(8,0.0);
  shared_ptr<matrixReal> T, U;
  tie(T,U) = R.svd(svd_vals);
  cout << *T << *U;
  T->setSub(1,0,*O.transpose());
  cout << *T;

  std::shared_ptr<matrixReal> H1(new matrixReal(2,2));
  std::shared_ptr<matrixReal> Sp1(new matrixReal(2,2));
  std::shared_ptr<matrixReal> Sz1(new matrixReal(2,2));
  Sp1->element(0,1) = 1.0;
  Sz1->element(0,0) = 0.5;
  Sz1->element(1,1) = -0.5;
  block testBlock(1,2,H1, Sp1, Sz1);
  Sz1->printMem();
  testBlock.enlarge(*H1, *Sp1, *Sz1);
  cout << "Hamiltonian:" << endl << *testBlock.H;
  cout << "Sp:" << endl << *testBlock.Sp;
  cout << "Sz:" << endl << *testBlock.Sz;
  Sz1->printMem();
  (*Sp1).~matrixReal();
  Sz1->printMem();

  return 0;
}
