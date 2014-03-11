//Given a Hamiltonian H, chain length L, vector set size m,
//computes the energy of the chain using the inifinite DMRG algorithm
//JL and JU are the lower and upper site coupling matrices
int dmrgInfiniteSystem(matrixReal &H, matrixReal &JU, matrixReal &JL, int L, int m);

class Block{
public:
	int nsites; //chain/number of sites
	int basisSize; //number of m's to keep
	std::vector<matrixReal> Ops; //operators for single site (initially)
}