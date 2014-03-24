#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
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
#if 1
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

#if 0
  //these are the initial matrices for H, Sp, and Sz, and WILL be changed by the enlarge function
  std::shared_ptr<matrixReal> H1(new matrixReal(2,2));
  auto sH1 = make_shared<matrixReal>(*H1);
  auto eH1 = make_shared<matrixReal>(*H1);

  std::shared_ptr<matrixReal> Sp1(new matrixReal(2,2));
  Sp1->element(0,1) = 1.0;
  auto sSp1 = make_shared<matrixReal>(*Sp1);
  auto eSp1 = make_shared<matrixReal>(*Sp1);

  std::shared_ptr<matrixReal> Sz1(new matrixReal(2,2));
  Sz1->element(0,0) = 0.5;
  Sz1->element(1,1) = -0.5;
  auto sSz1 = make_shared<matrixReal>(*Sz1);
  auto eSz1 = make_shared<matrixReal>(*Sz1);




  //Implemented in shoddy, procedural style. Will improve later, once the math all works.

  int maxKeepNum=20; //# maximum of eigenstates to keep
  int maxChainLen=10;

  block sysBlock(1, 2, maxKeepNum, sH1, sSp1, sSz1);
  block envBlock(1, 2, maxKeepNum, eH1, eSp1, eSz1);

  double error, energy;

  for (int currentLen = 0; currentLen < maxChainLen; currentLen++){
    cout << "System Length = " << 2*sysBlock.nSites+2 << endl;

    tie(error,energy)=sysBlock.infiniteDMRGStep(envBlock);

    cout << "Truncation Error : " << error << endl;
    cout << "E/L = " << std::setprecision(10) << energy/(sysBlock.nSites*2) << endl;
    sysBlock.H->printMem();
    cout << "******************************************************************************" << endl;
  }
#endif
  return 0;
}
