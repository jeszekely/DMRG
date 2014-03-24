#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "davidson.hpp"

using namespace std;

Davidson::Davidson(const int nV, const int size, const int guess, vector<double>& HD, function<matrixReal(matrixReal&)> h) : nVec(max(nV,2)), sizeVec(size), guessVec(guess), H(h), HDiag(HD){}

void Davidson::initialize(matrixReal& V)
{
  V.random();
  orthonorm(V);
  for (int ii = 0; ii < V.nc(); ii++) cout << V.dot(ii,ii) << endl;
  return;
}

void Davidson::orthonorm(matrixReal& V)
{
  //Orthogonalize via Gram-Schmidt
  for (int ii = 1; ii < V.nc(); ii++)
   {
    for (int jj = 0; jj < ii; jj++)
     {
      const double factor = V.dot(ii,jj);
        V.setSub(0,ii,*V.getSub(0,ii,V.nc(),1)-*V.getSub(0,jj,V.nc(),1)*factor);
     }
     const double norm =  V.dot(ii,ii);
     V.scaleCol(1.0/sqrt(norm),ii);
   }
  //Normalize
  for (int ii = 0; ii < V.nc(); ii++)
    V.scaleCol(1.0/sqrt(V.dot(ii,ii)),ii);
  return;
}

//Diagonalize matrix, return eigenvectors as columns of a matrixReal
std::shared_ptr<matrixReal> Davidson::diagonalize(vector<double>& outVals)
{
  int nevf = 0;
  vector<int> isConv(nVec,0);

//Initialize the guess matrix
  auto T = make_shared<matrixReal>(sizeVec,guessVec);
  initialize(*T);
  auto V = make_shared<matrixReal>(*T);
//Compute the subspace eigenpairs

  //vector <double> * PrevEigVals;
  for (int mm = guessVec; mm < maxIter; mm+=guessVec)
  {
    auto Hk = make_shared<matrixReal>(T->nc(),T->nc());
    auto Wk = make_shared<matrixReal>(sizeVec,T->nc());
    for (int ii = 0; ii < T->nc(); ii++)
    {
      Wk->setSub(0,ii,H(*V->getSub(0,ii,V->nr(),1)));
      for (int jj = 0; jj <= ii; jj++)
      {
        Hk->element(ii,jj) = Hk->element(jj,ii) = V->dot(*Wk->getSub(0,ii,V->nr(),1),jj,0);
      }
    }
    vector<double> eigVals(Hk->nc(),0.0);
    Hk->diagonalize(eigVals.data());
    matrixReal psi = (*V)*(*Hk);
    *Wk *= *Hk;

//Calculate the residuals, check for convergence
    auto R = make_shared<matrixReal>(*V);
    for(int ii = 0; ii<R->nc(); ii++)
    {
      matrixReal ri = *psi.getSub(0,ii,psi.nr(),1)*eigVals.at(ii) - *Wk->getSub(0,ii,Wk->nr(),1);
      R->setSub(0,ii,ri);
      if (sqrt(R->dot(ii,ii)) < tolerance)
      {
        T->setSub(0,ii,*psi.getSub(0,ii,T->nr(),1));
        outVals[ii] = eigVals[ii];
        isConv[ii] = 1;
      }
    }

//If all the requested vector are converged, exit for loop
    nevf = accumulate(isConv.begin(),isConv.begin()+isConv.size(),0);
    if (nevf >= nVec)
    {
      cout << "all is converged" << endl;
      break;
    }
//Calculate the correction vectors
    auto D = make_shared<matrixReal>(*R);
    for (int ii = 0; ii < D->nc(); ii++)
    {
      matrixReal Di(D->nr(),1);
      for (int jj = 0; jj < Di.nr(); jj++)
        Di(jj,0) = R->element(jj,ii)/(eigVals[ii]-HDiag[jj]);
      D->setSub(0,ii,Di);
    }

  //auto V = make_shared<matrixReal>(*T);
    auto Tnew = make_shared<matrixReal>(T->nr(), T->nc()+D->nc());
    Tnew->setSub(0,0,*V);
    Tnew->setSub(0,V->nc(),*D);
    *T = *Tnew;
    //Define D, append to set of vectors
  }


  return T;
}



