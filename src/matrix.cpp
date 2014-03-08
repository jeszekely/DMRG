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

#include "mkl.h"
#include "matrix.hpp"

using namespace std;

matrixReal::matrixReal(const int nr, const int nc) : matrixBase<double>(nr,nc){}
matrixReal::matrixReal(const matrixReal& o) : matrixBase<double>(o){}
matrixReal::matrixReal(matrixReal&& o ) : matrixBase<double>(std::move(o)){}

matrixReal& matrixReal::operator=(const matrixReal& o) {
  assert(nrows == o.nrows && ncols == o.ncols);
  copy_n(o.data(), o.size(), data());
}

matrixReal& matrixReal::operator*=(const matrixReal& o)
{
  *this = *this * o;
  return *this;
}

matrixReal matrixReal::operator*(const matrixReal& o) const
{
  const int m = nrows;
  const int k = ncols;
  const int n = o.ncols;
  assert(k == o.nrows);

  double alpha = 1.0;
  double beta = 0.0;

  matrixReal out(m, n);

  // it'd probably still be better to have an overloaded function that takes care of all the converting to memory locations stuff
  dgemm("N", "N", &m, &n, &k, &alpha, data(), &m, o.data(), &n, &beta, out.data(), &m);
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
      //          m, m, m, alpha, vals, m, o.vals, m, beta, vals, m);

  return out;
}


 	//  */*= (mkl);
    //  SVD/EigenvalueDecomp (mkl)
    //  isValid
    //  checkHermitian
    //  *sparsify
