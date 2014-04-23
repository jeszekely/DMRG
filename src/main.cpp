#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>

#include "matrix.hpp"
#include "block.hpp"
#include "dmrg.hpp"
#include "vector.hpp"
#include "davidson.hpp"

#include "input_parser.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
  programInputs Args("inputs.json");

//Test code for the matrix class
if (Args.TC_BasicMatrix_Run)
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
}

//test code for the davidson diagonalization
if (Args.TC_Davidson_Run)
{
  vectorMatrix Vecs(10,10);
  Vecs.random();
  cout << Vecs;
  for (int ii = 0; ii < Vecs.nc(); ii++)
    cout << Vecs.dot(ii,ii) << endl;
  Vecs.orthonormAll();
  cout << Vecs;
  for (int ii = 0; ii < Vecs.nc(); ii++)
    cout << Vecs.dot(ii,ii) << endl;

  vectorMatrix R(Args.TC_Davidson_MatrixSize,Args.TC_Davidson_MatrixSize);
  R.random();
  R += *R.transpose();
  for (int i = 0; i < Args.TC_Davidson_MatrixSize; ++i)
    R(i,i) += i*i*1.;

  auto RVecs = make_shared<vectorMatrix>(R);

  vector <double> vals(Args.TC_Davidson_MatrixSize, 0.0);

  clock_t t1, t2;
  double FullDiagTime, DavidsonTime;
  t1 = clock();
  R.diagonalize(vals.data());
  t2 = clock();
  FullDiagTime = (float(t2)-float(t1))/CLOCKS_PER_SEC;

  vector<double> RVals;
  t1 = clock();
  tie(RVecs,RVals) = diagonalizeDavidson(*RVecs,Args.TC_Davidson_nVecs, Args.TC_Davidson_nVecs, 100, Args.TC_Davidson_tolerance, &Args); //RDave.diagonalize();
  t2 = clock();
  DavidsonTime = (float(t2)-float(t1))/CLOCKS_PER_SEC;

  cout << "full eig:     " << setw(22) << setprecision(16) << vals[0] << endl;
  cout << "davidson eig: " << setw(22) << setprecision(16) << RVals[0] << endl;
  cout << endl;
  cout << "full eig:     " << setw(22) << setprecision(16) << vals[1] << endl;
  cout << "davidson eig: " << setw(22) << setprecision(16) << RVals[1] << endl;
  cout << "Code execution time:" << endl;
  cout << "full diag: " << FullDiagTime << " s" << endl << "davidson: " << DavidsonTime << " s" << endl;
}

//DMRG test code
if (Args.TC_DMRGSpin_Run)
{
  int maxChainLen=Args.TC_DMRGSpin_ChainLen; //This means the entire thing is 20, so 10 on each side
  vector<int> maxKeepNum = {Args.TC_DMRGSpin_TestTrunc1, Args.TC_DMRGSpin_TestTrunc2, Args.TC_DMRGSpin_TestTrunc3};
  finiteSystem finiteChain(maxChainLen, maxKeepNum);
  finiteChain.sweep();
}

  return 0;
}
