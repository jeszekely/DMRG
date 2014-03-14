//Class for calculating eigenvalues/eigenvectors using the Davidson algorithm
class Davidson
{
protected:
  const int nVec;
  const int sizeVec;
  std::function<matrixReal(matrixReal&)> H;
  std::vector<std::shared_ptr<matrixReal>> Vecs;

  //fill Vecs with proper starting guess
  void initialize();
public:
  //accepts the number of vectors, length of the vectors, and a lambda function describing the application of H
  Davidson(int nV, int size, std::function<matrixReal(matrixReal&)> h);

  //normalize the ith column of Vecs
  void normalize(int i);

  //Orthonormalize Vecs
  void orthonorm();

  std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>> diagonalize();
};

// Davidson algorithm:
// Choose u1 with ||u1|| = 1, U1 = [u1]
//for j = 1,2,...
//  w^j = Au^j
//  for k = 1,2,..j-1
//    bkj = (u^k)^H w^j
//    bjk = (u^j)^H w^k
//  bjj = (u^j)^H w^J
//  Compute largest eigenvalue of B, s, with eigenvector S, ||S|| = 1
//  y = UjS
//  y = Ay - sy
//  t = (DA - sI)^-1 r
//  t = t - UjUj^H t
//  u^(j+1) = t/||t||
//  U(j+1) = [Uj,u^(j+1)]