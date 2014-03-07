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

matrixReal::matrixReal(const int nr, const int nc) : matrixBase<double>(nr,nc){}
matrixReal::matrixReal(const matrixReal& o) : matrixBase<double>(o){}
matrixReal::matrixReal(matrixReal&& o ) : matrixBase<double>(std::move(o)){}

matrixReal& matrixReal::operator*=(const matrixReal& o)
{
	
	const int m=nrows;
	const int n=o.nrows;
	const int p=o.ncols;
	assert (m==n && n==p);

	double alpha=1;
	double beta=0;
	dgemm("N", "N",&m,&m,&m,&alpha,vals,&m,o.vals,&m,&beta,vals,&m);
    //cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
      //          m, m, m, alpha, vals, m, o.vals, m, beta, vals, m);

}


 	//  */*= (mkl);
    //  SVD/EigenvalueDecomp (mkl)
    //  isValid
    //  checkHermitian
    //  *sparsify