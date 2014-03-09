#include <iostream>
#include <cmath>
#include <complex>
#include <new>
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

// #include <boost/algorithm/string.hpp>
// #include <boost/property_tree/ptree.hpp>
// #include <boost/property_tree/json_parser.hpp>
// #include <boost/lexical_cast.hpp>

using namespace std;

typedef std::complex<double> cplx;

int main(int argc, char const *argv[])
{
	matrixReal M(5,6);
	M.makeIdentity();
	M(2,1) = 10.0;

	matrixReal N(6,2);
	N.makeIdentity();
	N(4,1) = -3.0;
	N(2,0) = -20.0;

	cout << M << N;

	matrixReal O = M*N;
	cout << O;

	matrixReal Q = N-N*2.0;
	cout << Q;
	return 0;
}
