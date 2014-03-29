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

// //Class for calculating eigenvalues/eigenvectors using the Davidson algorithm
// class Davidson
// {
// protected:
//   const int nVec;
//   const int sizeVec;
//   const int guessVec;
//   const double tolerance = 1.0e-8;
//   const int maxIter = 1000;
//   std::function<matrixReal(matrixReal&)> H;
//   std::vector<double> HDiag;
//   //std::vector<std::shared_ptr<matrixReal>> Vecs;

//   //fill Vecs with proper starting guess
//   void initialize(matrixReal& V);
// public:
//   //accepts the number of vectors, length of the vectors, and a lambda function describing the application of H
//   Davidson(const int nV, const int size, const int guess, std::vector<double>&, std::function<matrixReal(matrixReal&)> h);

//   //Orthonormalize Vecs
//   void orthonorm(matrixReal& V);

//   std::shared_ptr<matrixReal> diagonalize(std::vector<double> &);
// };
