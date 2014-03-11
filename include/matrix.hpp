#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>

template <typename T> class matrixBase
{
protected:
    size_t nrows, ncols;
    std::unique_ptr<T[]> vals;

public:
//  Constructor
  matrixBase(const int nr, const int nc) : nrows(nr), ncols(nc), vals(std::unique_ptr<T[]>(new T[nc*nr]))
  {
    zero();
  }

//  Copy Constructor
  matrixBase(const matrixBase& o) : nrows(o.nrows), ncols(o.ncols), vals(std::unique_ptr<T[]>(new T[nrows*ncols]))
  {
    std::copy_n(o.vals.get(), nrows*ncols, vals.get());
  }

//  Move Constructor
  matrixBase(matrixBase&& o) : nrows(o.nrows), ncols(o.ncols), vals(std::move(o.vals)) { };

//  Access functions
  size_t size() const { return nrows * ncols; }

  T* data() { return vals.get(); }
  const T* data() const { return vals.get(); }

//  Fill with zeroes
  void zero()
  {
    std::fill_n(vals.get(), nrows*ncols, T(0.0));
  }

//  Accessor functions
  T& element(const int row, const int col)
  {
    return vals[col*nrows + row];
  }

  const T& element(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

  T& operator()(const int row, const int col)
  {
    return element(row,col);
  }

  const T& operator()(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

//  Set only diagonal elements to unity
  void makeIdentity()
  {
    zero();
    for (int ii = 0; ii < std::min(ncols,nrows); ii++)
      vals[ii+ii*nrows] = T(1.0);
  }

//  Compute the trace
  T trace()
  {
    T sum = T(0.0);
    for (int ii = 0; ii < std::min(int(ncols),int(nrows)); ii++)
      sum += element(ii,ii);
    return sum;
  }

//  Apply a scaling factor to all elements
  void scale(const T a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
  }

//Extract a portion of the matrix starting at element(r,c), get (nr x nc matrix)
  matrixBase <T> getSub(int r, int c, int nr, int nc) const
  {
    assert(r + nr < nrows && c + nc < ncols);
    matrixBase <T> out(nr,nc);
    for (int ii = 0; ii < nr; ii++)
    {
      for (int jj = 0; jj < nc; jj++)
      {
        out(ii,jj) = element(ii+r,jj+c);
      }
    }
    return out;
  }

//Place matrix o at position (r,c)
  template <typename U>
  void setSub(int r, int c, U o)
  {
    assert(r + o.nrows <= nrows && c + o.ncols <= ncols);
    for (int ii = 0; ii < o.nrows; ii++)
    {
      for (int jj = 0; jj < o.ncols; jj++)
      {
        element(ii+r,jj+c) = o(ii,jj);
      }
    }
    return;
  }

  template <typename U> friend std::ostream &operator<<(std::ostream &out, const matrixBase <U> &o);
};

class matrixReal : public matrixBase<double>
{
public:
  matrixReal(const int nr, const int nc);
  matrixReal(const matrixReal&);
  matrixReal(matrixReal&&);

//Matrix-Matrix operations
  matrixReal& operator=(const matrixReal&);
  matrixReal operator*(const matrixReal&) const;
  matrixReal& operator*=(const matrixReal&);
  matrixReal operator+(const matrixReal&) const;
  matrixReal& operator+=(const matrixReal&);
  matrixReal operator-(const matrixReal&) const;
  matrixReal& operator-=(const matrixReal&);

//Scalar-Matrix Operations
//Note: binary scalar operations only work as rhs operators at the moment
  matrixReal operator*(const double&) const;
  matrixReal operator/(const double&) const;
  matrixReal& operator*=(const double&);
  matrixReal& operator/=(const double&);

//BLAS and LAPACK routines
//Diagonalize matrix, place eigenvalues in the vector prodived
//NOTE: Assumes a symmetric matrix
  void diagonalize(double* eigVals);
  std::shared_ptr<matrixReal> transpose() const;
  std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>>svd(std::vector<double>&);

//Compute the kronecker product of two matrices
//Note: This will be used as a shortcut to create a superblock matrix for small site matrices
//This should NOT be used for larger systems as it is very memory intensive
  matrixReal kron(matrixReal &o) const;


};

//Overload the << operator to print a matrix
template <typename T>
std::ostream &operator<<(std::ostream &out, const matrixBase <T> &o)
{
  for (int row = 0; row < std::min(10,int(o.nrows)); row++)
  {
    for (int col = 0; col < std::min(10,int(o.ncols)); col++)
    {
      out << std::setprecision(3) << o(row,col) << "\t";
    }
    out << "\n";
  }
  out << std::endl;
  return out;
};

//class matrixComplex

