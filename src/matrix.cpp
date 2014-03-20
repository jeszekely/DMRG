#include <iostream>
#include <cmath>
#include <complex>
#include <memory>
#include <algorithm>
#include <stdexcept>

#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"
#include "utilities.hpp"

using namespace std;

matrixReal::matrixReal(const int nr, const int nc) : matrixBase<double>(nr,nc){}
matrixReal::matrixReal(const matrixReal& o) : matrixBase<double>(o){}
matrixReal::matrixReal(matrixReal&& o ) : matrixBase<double>(std::move(o)){}

matrixReal& matrixReal::operator=(const matrixReal& o)
{
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
  dgemm_("N","N", nrows, o.ncols, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), nrows);
  return out;
}

//vector dot product
double matrixReal::operator%(const matrixReal& o) const
{
  assert(ncols == 1 && o.ncols == 1 && nrows == o.nrows);
  return ddot_(nrows, data(), 1, o.data(), 1);
}

//vector dot product where vectors are two columns in matrix
double matrixReal::dot(const int ii, const int jj) const
{
  assert(ncols > max(ii,jj));
  return ddot_(nrows, &element(0,ii), 1, &element(0,jj), 1);
}

//vector dot product where vectors are two columns in separate matricies
double matrixReal::dot(const matrixReal& o, const int ii, const int jj) const
{
  assert(ncols > ii && o.ncols > jj && nrows == o.nrows);
  return ddot_(nrows, &element(0,ii), 1, &o.element(0,jj), 1);
}

//Matrix multiplication, left matrix is transposed
matrixReal matrixReal::operator|(const matrixReal& o) const
{
  assert(nrows == o.nrows);
  matrixReal out(ncols, o.ncols);
  dgemm_("T","N", ncols, o.ncols, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), ncols);
  return out;
}

//Matrix multiplication, right matrix is transposed
matrixReal matrixReal::operator^(const matrixReal& o) const
{
  assert(ncols == o.ncols);
  matrixReal out(nrows, o.nrows);
  dgemm_("N","T", nrows, o.nrows, o.nrows, 1.0, data(), nrows, o.data(), o.nrows, 0.0, out.data(), nrows);
  return out;
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
  *this = *this + o;
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
  *this = *this - o;
  return *this;
}

matrixReal matrixReal::operator*(const double& a) const
{
  matrixReal out(*this);
  out *= a;
  return out;
}

matrixReal& matrixReal::operator*=(const double& a)
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

matrixReal& matrixReal::operator/=(const double& a)
{
  scale(1.0/a);
  return *this;
}

void matrixReal::diagonalize(double* eigVals)
{
  assert (nrows == ncols);
  int info;
  int lwork = -1;
  double wkopt;
  dsyev_("V", "U", nrows, data(), nrows, eigVals, &wkopt, lwork, info);
  lwork = int(wkopt);
  std::unique_ptr <double[]> work (new double [lwork]);
  dsyev_("V", "U", nrows, data(), nrows, eigVals, work.get(), lwork, info);
  if (info > 0)
    throw std::runtime_error("Unable to diagonalize matrix");
  return;
}

std::shared_ptr<matrixReal> matrixReal::transpose() const
{
  auto out = make_shared<matrixReal> (ncols,nrows);
  mkl_domatcopy_("C","T",nrows,ncols,1.0,data(),nrows,out->data(),ncols);
  return out;
}

tuple<shared_ptr<matrixReal>, shared_ptr<matrixReal>> matrixReal::svd(vector<double>& s)
{
  assert(s.size() >= std::min(ncols,nrows));
  auto u = make_shared<matrixReal>(nrows, nrows);
  auto vT = make_shared<matrixReal>(ncols, ncols);

  int lwork = -1;
  int info = 0;
  dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, s.data(), lwork, info);
  lwork = s[0];
  if (lwork <= 0)
      throw runtime_error("dgesvd faied allocating lwork value");
  unique_ptr<double[]> work(new double[lwork]);
  dgesvd_("A","A", nrows, ncols, data(), nrows, s.data(), u->data(), nrows, vT->data(), ncols, work.get(), lwork, info);
  if (info != 0)
      throw runtime_error("dgesvd faied matrix decomposition");
  return make_tuple(u,vT);
}

matrixReal matrixReal::kron(matrixReal &o) const
{
  matrixReal out(nrows*o.nrows,ncols*o.ncols);
  for (int ii = 0; ii < nrows; ii++)
  {
    for (int jj = 0; jj < ncols; jj++)
    {
      out.setSub(ii*o.ncols,jj*o.nrows,o*element(ii,jj));
    }
  }
  return out;
}

  //  Hermitian Conjugate
  //  *sparsify
