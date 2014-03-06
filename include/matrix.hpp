#include <memory>

template <typename T> class matrix
{
protected:
	int nrows, ncols;
	std::unique_ptr<T[]> vals;

public:
  // 	Constructor
	matrix(int nr, int nc) : nrows(nr), ncols(nc), vals(std::unique_ptr<T[]>(new T[nc*nr]))
	{
    zero();
	};

  matrix(const matrix& o) : nrows(o.nrows), ncols(o.ncols), vals(std::unique_ptr<T[]>(new T[nrows*ncols])) {
    std::copy_n(o.vals.get(), nrows*ncols, vals.get());
  }

  matrix(matrix&& o) : nrows(o.nrows), ncols(o.ncols), vals(std::move(o.vals)) { }

  void zero() {
    std::fill_n(vals.get(), nrows*ncols, T(0.0));
  }

  T& element(const int row, const int col) {
    return vals[col + row*ncols];
  }

  const T& element(const int row, const int col) const {
    return vals[col + row*ncols];
  }

  T& operator()(const int row, const int col) {
    return element(row,col);
  }

  const T& operator()(const int row, const int col) const {
    return vals[col + row*ncols];
  }

  //	get a value
	T get(int row, int col)
	{
		return vals[col + row*ncols];
	};

  //	set a particular value
	int set(int row, int col, T val)
	{
		vals[col+row*ncols] = val;
		return 0;
	};

  //	Set only diagonal elements to unity
	void makeIdentity()
	{
		for (int ii = 0; ii < std::min(ncols,nrows); ii++)
		{
			vals[ii+ii*ncols] = (T)1.0;
		}
	};

	void printMatrix()
	{
		for (int col = 0; col < std::min(10,nrows); col++)
		{
			for (int row = 0; row < std::min(10,ncols); row++)
			{
//				cout << this->get(row,col) << " ";
				std::cout << "(" << row << " " << col << " " << get(row,col) << ")" << " ";
			}
			std::cout << std::endl;
		}
	};

	// 	Move
	// 	Copy
	// 	*/*= (mkl);
	// 	+/+=
	// 	SVD/EigenvalueDecomp (mkl)
	// 	Trace
	// 	isValid
	// 	checkHermitian
	// 	*sparsify


	/* data */
};


