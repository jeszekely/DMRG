//Given a Hamiltonian H, chain length L, vector set size m,
//computes the energy of the chain using the inifinite DMRG algorithm
//JL and JU are the lower and upper site coupling matrices
class block{
public:
	int nSites; //chain/number of sites
	int basisSize; //number of m's to keep
	std::shared_ptr<matrixReal> H;
	std::shared_ptr<matrixReal> JU;
	std::shared_ptr<matrixReal> JL;
	std::vector<matrixReal *> Ops; //operators for single site (initially)

	block(int, int, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>);
	void enlarge(block &o);
};

int dmrgInfiniteSystem(block& system, int L, int m);
