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

 	//  */*= (mkl);
    //  SVD/EigenvalueDecomp (mkl)
    //  isValid
    //  checkHermitian
    //  *sparsify