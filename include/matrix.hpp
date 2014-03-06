using namespace std; 

template <typename T> class matrix
{
public:
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
		for (int ii = 0; ii < min(ncols,nrows); ii++)
		{
			vals[ii+ii*ncols] = (T)1.0; 
		}
	};

	void printMatrix()
	{
		for (int col = 0; col < min(10,nrows); col++)
		{
			for (int row = 0; row < min(10,ncols); row++)
			{
//				cout << this->get(row,col) << " "; 
				cout << "(" << row << " " << col << " " << this->get(row,col) << ")" << " "; 
			}
			cout << endl;
		}
	};

// 	Destructor
	~matrix()
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


	/* data */
};


