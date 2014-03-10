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
}

matrixReal matrixReal::transpose() const
{
  matrixReal out(ncols,nrows);
  mkl_domatcopy_("C","T",nrows,ncols,1.0,data(),nrows,out.data(),ncols);
  return out;
}



matrixReal matrixReal::svd()
{
  //Uses LAPACKE_dgesvd to get the singular value decomposition
  //Right now returns 'data', which is modified by this function, 
  //which is probably not what we want. Fix this.
  

  double * superb = new double [min(nrows,ncols)-1];  
  lapack_int m = nrows, n = ncols;

  matrixReal s(ncols,1);
  matrixReal u(nrows,nrows);
  matrixReal vt(ncols,ncols);

  lapack_int ldu = nrows, ldvt = ncols, info, lwork, lda=nrows;  


  info = LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'A', 'A', nrows, ncols, data(), lda,  
               s.data(), u.data(), ldu, vt.data(), ldvt, superb );
 
 if( info > 0 ) {
    printf( "The algorithm computing SVD failed to converge.\n" );
    exit( 1 );
  }  

  //Can print out the SVD results if you want
  cout << "Singular Values" << endl << s;
  cout << "Left singular vectors (stored columnwise)" << endl << u;
  cout << "Right singular vectors (stored rowwise)" << endl << vt;

  return *this;
}

  //  Invert
  //  kron product
  //  Remove a row or column
  //  Extract submatrix
  //  Extract one row or column (derived from above)
  //  Transpose / Hermitian Conjugate
  //  SVD/EigenvalueDecomp (mkl)
  //  separate eigenvalue function?
  //  isValid
  //  *sparsify
