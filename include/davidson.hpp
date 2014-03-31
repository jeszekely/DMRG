#ifndef DMRG_DAVIDSON
#define DMRG_DAVIDSON
#include "matrix.hpp"
#include "vector.hpp"

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
public:
  Davidson(genMatrix H, int initialVecs, int numVecs, int maxIterations, double error);
  std::tuple<std::shared_ptr<vectorMatrix>,std::vector<double>> diagonalize();
};

#endif