#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"
#include "davidson.hpp"

using namespace std;

Davidson::Davidson(const int nV, const int size, function<matrixReal(matrixReal&)> h) : nVec(min(2,nV)), sizeVec(size), H(h), Vecs(vector<shared_ptr<matrixReal>>(nV))
{
  for (shared_ptr<matrixReal>& v : Vecs)
    v = make_shared<matrixReal>(sizeVec,1);
  initialize();
}

void Davidson::initialize()
{
  for (shared_ptr<matrixReal>& v : Vecs)
  {
    //Fill with random values then normalize to 1.0
    //This is probably not ideal for a starting guess to the Davidson algorithm
    v->random();
    v->scale(1.0/sqrt((*v|*v)(0,0)));
  }
  orthonorm();
  for (shared_ptr<matrixReal>& v : Vecs) cout << (*v|*v)(0,0) << endl;
  return;
}

void Davidson::orthonorm()
{
  //Orthogonalize via Gram-Schmidt
  for (int ii = 1; ii < nVec; ii++)
   {
    for (int jj = 0; jj < ii; jj++)
     {
      const double factor = (*Vecs.at(ii)|*Vecs.at(jj))(0,0);
       *(Vecs.at(ii)) -= *Vecs.at(jj)*factor;
     }
     const double norm =  (*Vecs.at(ii)|*Vecs.at(ii))(0,0);
     Vecs.at(ii)->scale(1.0/sqrt(norm));
   }
  //Normalize
  for (shared_ptr<matrixReal>& v : Vecs)
    v->scale(1.0/sqrt((*v|*v)(0,0)));
  return;
}

//Diagonalize matrix, return eigenvectors as columns of a matrixReal
std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>> Davidson::diagonalize()
{
  for (shared_ptr<matrixReal>& v : Vecs)
  {
      cout << *v << endl;
  }
}


// Davidson algorithm:
// Choose u1 with ||u1|| = 1, U1 = [u1]
//for j = 1,2,...
//  w^j = Au^j
//  for k = 1,2,..j-1
//    bkj = (u^k)^H w^j
//    bjk = (u^j)^H w^k
//  bjj = (u^j)^H w^J
//  Compute largest eigenvalue of B, s, with eigenvector S, ||S|| = 1
//  y = UjS
//  y = Ay - sy
//  t = (DA - sI)^-1 r
//  t = t - UjUj^H t
//  u^(j+1) = t/||t||
//  U(j+1) = [Uj,u^(j+1)]


