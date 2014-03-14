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

Davidson::Davidson(const int nV, const int size, function<matrixReal(matrixReal&)> h) : nVec(nV), sizeVec(size), H(h), Vecs(vector<shared_ptr<matrixReal>>(nV))
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
       *(Vecs.at(ii)) -= (*(Vecs.at(jj)))*((*Vecs.at(ii)|*Vecs.at(jj))(0,0) / (*Vecs.at(jj)|*Vecs.at(jj))(0,0));
     }
   }
  //Normalize
  for (shared_ptr<matrixReal>& v : Vecs)
    v->scale(1.0/sqrt((*v|*v)(0,0)));
  return;
}


