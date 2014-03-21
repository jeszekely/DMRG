//Class for calculating eigenvalues/eigenvectors using the Davidson algorithm
class Davidson
{
protected:
  const int nVec;
  const int sizeVec;
  const int guessVec;
  const double tolerance = 1.0e-8;
  const int maxIter = 1000;
  std::function<matrixReal(matrixReal&)> H;
  std::vector<double> HDiag;
  //std::vector<std::shared_ptr<matrixReal>> Vecs;

  //fill Vecs with proper starting guess
  void initialize(matrixReal& V);
public:
  //accepts the number of vectors, length of the vectors, and a lambda function describing the application of H
  Davidson(const int nV, const int size, const int guess, std::vector<double>&, std::function<matrixReal(matrixReal&)> h);

  //Orthonormalize Vecs
  void orthonorm(matrixReal& V);

  std::shared_ptr<matrixReal> diagonalize(std::vector<double> &);
};
