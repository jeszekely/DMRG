#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "davidson.hpp"

using namespace std;

genMatrix::genMatrix(size_t nr, size_t nc, std::function<vectorMatrix(vectorMatrix&)> h, vector<double> &Diag) : nrows(nr), ncols(nc), H(h),diags(Diag) {}

double genMatrix::diagElem(int ii) {return diags[ii];}

int genMatrix::nc() {return ncols;}

int genMatrix::nr() {return nrows;}

vectorMatrix genMatrix::operator*(vectorMatrix& o){return H(o);}

Davidson::Davidson(genMatrix h, int iVecs, int nVecs, int maxIts, double err) : H(h), initialVecs(iVecs), numVecs(nVecs), maxIterations(maxIts), tolerance(err) {}

tuple<std::shared_ptr<vectorMatrix>,vector<double>> Davidson::diagonalize()
{
  int n = H.nr(); //Size of a vector
  int mmax = n/2; //Maximum iterations

  //initial guess set to unit vectors
  vectorMatrix t(n,initialVecs);
  t.makeIdentity();

  //Array to hold guess vectors
  vectorMatrix V(n,n);

  //eigenvector and eigenvalue return arrays
  auto eigVecs = make_shared<vectorMatrix>(n,numVecs);
  vector<double> eigVals(numVecs,1.0);
  vector<double> theta_old(eigVals);

  vectorMatrix Ht = H*t; //matrix acting on test vectors

  for (int m = initialVecs; m < mmax; m+=initialVecs)
  {
    if (m <= initialVecs)
    {
      V.setSub(0,0,t);
    }
    else if (m > initialVecs) copy_n(eigVals.data(),numVecs,theta_old.data());
    vectorMatrix T = (t|Ht);
    vector<double> TeigVals(T.nc(),0.0);
    T.diagonalize(TeigVals.data());
    for (int jj = 0; jj<initialVecs; jj++)
    {
      vectorMatrix r = Ht*T - *V.getSub(0,0,V.nr(),m+1)*T*TeigVals[jj];
      eigVecs->setSub(0,jj,*V.getSub(0,0,V.nr(),m+1)*T);
      vectorMatrix q(r);
      for (int kk = 0; kk < q.nr(); kk++)
        q(kk,0) = (1.0/(TeigVals[jj] - H.diagElem(kk)))*q(kk,0);
      V.setSub(0,m+1+jj,q);
    }
    double norm = 0.0;
    for (int jj = 0; jj < numVecs; jj++)
      norm += abs(pow(theta_old[jj] - eigVals[jj],2));
    norm = sqrt(norm);
    if (norm < tolerance) break;
  }

  return make_tuple(eigVecs,eigVals);
}
