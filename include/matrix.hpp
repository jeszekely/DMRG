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

//  Fill with zeroes
    void zero()
    {
        std::fill_n(vals.get(), nrows*ncols, T(0.0));
    }

//  Accessor functions
    T& element(const int row, const int col) 
    {
        return vals[col + row*ncols];
    }

    const T& element(const int row, const int col) const 
    {
        return vals[col + row*ncols];
    }

    T& operator()(const int row, const int col) 
    {
        return element(row,col);
    }

    const T& operator()(const int row, const int col) const 
    {
        return vals[col + row*ncols];
    }

//  Set only diagonal elements to unity
    void makeIdentity()
    {
        zero();
        for (int ii = 0; ii < std::min(ncols,nrows); ii++)
            vals[ii+ii*ncols] = T(1.0);
    }

//  Print a small portion for error checking
    void printMatrix() const
    {
        for (int row = 0; row < std::min(10,int(nrows)); row++)
        {
            for (int col = 0; col < std::min(10,int(ncols)); col++)
            {
              std::cout << std::setprecision(2) << element(row,col) << " ";
            }
            std::cout << std::endl;
        }
    }

//  Compute the trace
    T trace()
    {
        T sum = T(0.0);
        for (int ii = 0; ii < std::min(int(ncols),int(nrows)); ii++)
            sum += element(ii,ii);
        return sum;
    }

    matrixBase& operator+=(const matrixBase& o)
    {
        assert(nrows == o.nrows);
        assert(ncols == o.ncols);
        for (int ii = 0; ii < ncols*nrows; ii++)
            vals[ii] += o.vals[ii];
        return *this; 
    }

    matrixBase& operator-=(const matrixBase& o)
    {
        assert(nrows == o.nrows);
        assert(ncols == o.ncols);
        for (int ii = 0; ii < ncols*nrows; ii++)
            vals[ii] -= o.vals[ii];
        return *this; 
    }
    //  +, -/-=

};

class matrixReal : public matrixBase<double>
{
public:
    matrixReal(const int nr, const int nc);
    matrixReal(const matrixReal&);
    matrixReal(matrixReal&&);
};

//class matrixComplex

