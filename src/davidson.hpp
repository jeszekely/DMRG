#ifndef DMRG_DAVIDSON
#define DMRG_DAVIDSON
#include "matrix.hpp"
#include "vector.hpp"
#include "input_parser.hpp"

//generalized matrix class designed to be used as input to the davidson algorithm
//contains the number of rows and cols in a matrix, along with a multiplication procedure
class genMatrix
{
protected:
  size_t nrows;
  size_t ncols;
  std::function<vectorMatrix(vectorMatrix&)> H;
  std::vector<double> diags;

public:
  vectorMatrix operator*(vectorMatrix&);
  genMatrix(size_t nr, size_t nc, std::function<vectorMatrix(vectorMatrix&)> H, std::vector<double>&);
  int nc();
  int nr();
  double diagElem(int ii);
};

class Davidson
{
protected:
  genMatrix H;
  int initialVecs;
  int numVecs;
  int maxIterations;
  double tolerance;
  bool verbose = false;
public:
  Davidson(genMatrix H, int initialVecs, int numVecs, int maxIterations, double error);
  std::tuple<std::shared_ptr<vectorMatrix>,std::vector<double>> diagonalize();
  void init(programInputs &P);
};

//Type T must be "matrix-like":
//  must have nr(), nc() functions
//  access to each element
//  multiplication operator
template <typename T>
std::tuple<std::shared_ptr<vectorMatrix>,std::vector<double>> diagonalizeDavidson(T &o, int iGuess, int nV, int max, double tol, programInputs *A = NULL)
{
  assert(o.nr() == o.nc());
  std::vector<double> eigVals;
  std::vector<double> oDiags(o.nr());
  for (int ii = 0; ii < o.nr(); ii++)
    oDiags[ii] = o(ii,ii);
  genMatrix oGen(o.nr(),o.nc(),[&o](vectorMatrix &vec){return o*vec;},oDiags);
  Davidson oDave(oGen, iGuess, nV, max, tol);
  if (A) oDave.init(*A);
  return oDave.diagonalize();
}


#endif