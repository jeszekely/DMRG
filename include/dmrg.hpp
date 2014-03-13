//Given a Hamiltonian H, chain length L, vector set size m,
//computes the energy of the chain using the inifinite DMRG algorithm
//JL and JU are the lower and upper site coupling matrices
class block{
public:
	int nSites; //chain/number of sites
	int basisSize; //number of m's to keep
	std::shared_ptr<matrixReal> H;
	std::shared_ptr<matrixReal> Sp;
	std::shared_ptr<matrixReal> Sz;
	std::vector<matrixReal *> Ops; //operators for single site (initially)

	block(int, int, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>);
	void enlarge(matrixReal &H1, matrixReal &Sp1, matrixReal &Sz1);
	void rotateTruncate(matrixReal& transformationMatrix, int maxSize);
};

std::shared_ptr<matrixReal> buildSuperblock(block& sysBlock, block& envBlock);

std::shared_ptr<matrixReal> makeReducedDM(matrixReal& groundWfxn);

std::shared_ptr<matrixReal> makeTransformationMatrix(matrixReal& reducedDM, int basisSize, int keepNum);

double truncationError(std::vector<double>& eigenvals, int basisSize, int keepNum);


int dmrgInfiniteSystem(block& system, int L, int m);
