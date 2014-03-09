#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>

#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <memory>

#include "mkl.h"
#include "matrix.hpp"

using namespace std;

matrixReal::matrixReal(const int nr, const int nc) : matrixBase<double>(nr,nc){}
matrixReal::matrixReal(const matrixReal& o) : matrixBase<double>(o){}
matrixReal::matrixReal(matrixReal&& o ) : matrixBase<double>(std::move(o)){}

matrixReal& matrixReal::operator=(const matrixReal& o) {
  assert(nrows == o.nrows && ncols == o.ncols);
  copy_n(o.data(), o.size(), data());
  return *this;
}

matrixReal& matrixReal::operator*=(const matrixReal& o)
{
  *this = *this * o;
  return *this;
}

matrixReal matrixReal::operator*(const matrixReal& o) const
{
  assert(ncols == o.nrows);
  matrixReal out(nrows, o.ncols);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
               nrows, o.ncols, o.nrows, 1.0, data(), o.nrows, o.data(), o.ncols, 0.0, out.data(), o.ncols);
  return out;

  // const int m = nrows;
  // const int k = ncols;
  // const int n = o.ncols;
  // double alpha = 1.0;
  // double beta = 0.0;
  // it'd probably still be better to have an overloaded function that takes care of all the converting to memory locations stuff
  //dgemm("N", "N", &m, &n, &k, &alpha, data(), &m, o.data(), &n, &beta, out.data(), &m);

}

matrixReal matrixReal:: operator+(const matrixReal& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixReal out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), plus<double>()); 
  return out;
}

matrixReal& matrixReal::operator+=(const matrixReal& o)
{
  *this = *this * o;
  return *this;
}

matrixReal matrixReal::operator-(const matrixReal& o) const
{
  assert(ncols == o.ncols && nrows == o.nrows);
  matrixReal out(nrows, o.ncols);
  transform(data(), data()+size(), o.data(), out.data(), minus<double>()); 
  return out;

}

matrixReal& matrixReal::operator-=(const matrixReal& o)
{
  *this = *this * o;
  return *this;
}

matrixReal matrixReal::operator*(const double& a) const
{
  matrixReal out(*this);
  out *= a;
  return out;
}


matrixReal& matrixReal::operator *= (const double& a)
{
  scale(a);
  return *this;
}

matrixReal matrixReal::operator/(const double& a) const
{
  matrixReal out(*this);
  out /= a;
  return out;
}

matrixReal& matrixReal::operator /= (const double& a)
{
  scale(1.0/a);
  return *this;
}

    //  Determinant
    //  Invert
    //  Transpose / Hermitian Conjugate
    //  SVD/EigenvalueDecomp (mkl)
    //    separate eigenvalue function?
    //  isValid
    //  checkHermitian
    //  *sparsify
