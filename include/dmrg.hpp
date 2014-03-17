//Given a Hamiltonian H, chain length L, vector set size m,
//computes the energy of the chain using the inifinite DMRG algorithm
//JL and JU are the lower and upper site coupling matrices
class block{
public:
	int nSites; //chain/number of sites
	int basisSize; //number of m's to keep
	int maxKeepNum; //maximum number of eigenvectors to keep
	std::shared_ptr<matrixReal> H;
	std::shared_ptr<matrixReal> Sp;
	std::shared_ptr<matrixReal> Sz;
	std::vector<matrixReal *> Ops; //operators for single site (initially)

	block(int, int, int, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>);
	
	void enlarge(matrixReal &H1, matrixReal &Sp1, matrixReal &Sz1);
	void rotateTruncate(matrixReal& transformationMatrix);
	void makeTransformationMatrix(matrixReal& reducedDM, block& environmentBlock);

	//Returns the ground state energy and the truncation error
	std::tuple<double, double> infiniteDMRGStep(block& environmentBlock);

	//Returns the ground state energy and wavefunction 
	std::tuple<double, std::shared_ptr<matrixReal>> buildSuperblock(block& envBlock);

	//Returns the truncation error and the diagonalized reduced density matrix
	std::tuple<double, std::shared_ptr<matrixReal>> makeReducedDM(matrixReal& groundWfxn);

};




int dmrgInfiniteSystem(block& system, int L, int m);
