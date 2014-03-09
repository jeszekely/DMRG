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

//  Compute the trace
    T trace()
    {
        T sum = T(0.0);
        for (int ii = 0; ii < std::min(int(ncols),int(nrows)); ii++)
            sum += element(ii,ii);
        return sum;
    }

    void scale(const T a)
    {
      std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
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
    matrixReal& operator *= (const double&);
    matrixReal& operator /= (const double&);

};

template <typename T>
std::ostream &operator<<(std::ostream &out, const matrixBase <T> &o)
{
    for (int row = 0; row < std::min(10,int(o.nrows)); row++)
    {
        for (int col = 0; col < std::min(10,int(o.ncols)); col++)
        {
          out << std::setprecision(2) << o(row,col) << " ";
        }
        out << "\n";
    }
    out << "\n";
    return out;
};

//class matrixComplex

