#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "davidson.hpp"
#include "vector.hpp"

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
  //t.makeIdentity();
  t.random();
  t.orthonormAll();

  //Array to hold guess vectors
  vectorMatrix V(n,n);
  //eigenvector and eigenvalue return arrays
  auto eigVecs = make_shared<vectorMatrix>(n,numVecs);
  vector<double> eigVals(numVecs,1.0);
  vector<double> theta_old(numVecs,100.0);
  vectorMatrix Ht(n,mmax);
  Ht.setSub(0,0,H*t); //matrix acting on test vectors

  for (int m = initialVecs; m < mmax; m+=initialVecs)
  {
    cout << "Step #" << m << endl;
    if (m <= initialVecs) {V.setSub(0,0,t);}
    else copy_n(eigVals.data(),numVecs,theta_old.data());

    vectorMatrix HtTemp = vectorMatrix(*V.getSub(0,0,V.nr(),m));
    Ht.setSub(0,0,H*HtTemp);
    vectorMatrix T = (*V.getSub(0,0,V.nr(),m)|HtTemp);
    T.diagonalize(eigVals.data());
    for (int jj = 0; jj<initialVecs; jj++)
    {
      vectorMatrix r = HtTemp*(*T.getSub(0,jj,T.nr(),1)) - *V.getSub(0,0,V.nr(),T.nr())*(*T.getSub(0,jj,T.nr(),1))*eigVals[jj];
      auto mat = make_shared<vectorMatrix>(*V.getSub(0,0,V.nr(),T.nr())*(*T.getSub(0,jj,T.nr(),1)));
      if (jj < numVecs) eigVecs->setSub(0,jj,*mat);
      vectorMatrix q(r);
      for (int kk = 0; kk < q.nr(); kk++)
        q(kk,0) = (1.0/(eigVals[jj] - H.diagElem(kk)))*q(kk,0);
      V.setSub(0,m+jj,q);
    }
    V.orthonormAll();
    double norm = 0.0;
    for (int jj = 0; jj < numVecs; jj++){
      norm += abs(pow(theta_old[jj] - eigVals[jj],2));}
    norm = sqrt(norm);
    if (norm < tolerance) break;
  }
  return make_tuple(eigVecs,eigVals);
}
