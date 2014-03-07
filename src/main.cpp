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
	matrix <float> M(5,6); 
	M.printMatrix();
	M.makeIdentity(); 
	std::cout << endl; 
	M.printMatrix();
	std::cout << endl; 
	M(2,1) = 10.0; 
	M.printMatrix();
	std::cout << "The trace is " << M.trace() << endl; 
	return 0;
}
