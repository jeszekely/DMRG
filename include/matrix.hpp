using namespace std; 

template <typename T> class matrix
{
private:
	int nrows, ncols;
	T *vals;
	
// 	Constructor
	matrix(int nr, int nc)
	{
		nrows 	= nr; 
		ncols 	= nc; 
		vals 	= new T [nc*nr]; 
		
//	zero all elements initially
		for (int ii = 0; ii < ncols*nrows; ii++)
		{
			vals[ii] = (T)0.0; 
		}
	};

//	get a value
	T get(int ii, int jj)
	{
		return vals[jj + ii*ncols]; 
	};
		
//	set a particular value
	int set(int ii, int jj, T val)
	{
		vals[jj+ii*ncols] = val; 
		return 0; 
	};

//	Set only diagonal elements to unity
	void makeIdentity()
	{
		for (int ii = 0; ii < min(ncols,nrows); ii++)
		{
			vals[ii+ii*ncols] = (T)1.0; 
		}
	};

	void printMatrix()
	{
		for (int ii = 0; ii < min(6,ncols); ii++)
		{
			for (int jj = 0; jj < min(6,nrows); jj++)
			{
				cout << this->get(ii,jj) << " "; 
			}
			cout << endl;
		}
	}

// 	Destructor
	void ~matrix()
	{
		delete [] vals; 
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

protected:

	matrix(arguments);
	~matrix();


	/* data */
};


