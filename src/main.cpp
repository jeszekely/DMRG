#include <iostream>
#include <cmath>
#include <complex>
#include <new>

#include <fstream>
#include <iomanip>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>
#include <memory>
#include <string>

#include "mkl.h"
#include "matrix.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

typedef std::complex<double> cplx;

int main(int argc, char const *argv[])
{
	matrix <double> M(5,2); 
	M.printMatrix();
	M.makeIdentity(); 
	M.printMatrix();
	return 0;
}
