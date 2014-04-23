//Given a Hamiltonian H, chain length L, vector set size m,
//computes the energy of the chain using the inifinite DMRG algorithm
//JL and JU are the lower and upper site coupling matrices


class infiniteSystem{
public:
	int sysLength;
	std::vector<int> sweepSizes; //for the infinite case it will only take the first value for this
	//Constructor
	infiniteSystem(int, std::vector<int>);

	void runInfinite();

	//Returns the ground state energy and the truncation error
	std::tuple<double, double> infiniteDMRGStep(block& sysBlock, block& envBlock, int maxEigenStates);

	std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>> makeSingleSiteOperators();

	void enlarge(block& b);

	//Returns the ground state energy and wavefunction
	std::tuple<double, std::shared_ptr<matrixReal>> buildSuperblock(block& sysBlock, block& envBlock);

	//Returns the truncation error and the diagonalized reduced density matrix
	std::tuple<double, std::shared_ptr<matrixReal>> makeReducedDM(int sysBasisSize, int envBasisSize, matrixReal& groundWfxn, int maxEigenStates);

	//void rotateTruncate(matrixReal& transformationMatrix);
	void makeTransformationMatrix(block& sysBlock, matrixReal& reducedDM, int maxEigenStates);

	void graphics(int sysBlockLength, int envBlockLength, bool direction);

};


class finiteSystem : public infiniteSystem{
public:

	//Constructor
	finiteSystem(int, std::vector<int>);

	std::vector<std::shared_ptr<block>> leftBlocks;
	std::vector<std::shared_ptr<block>> rightBlocks;

	std::tuple<std::vector<std::shared_ptr<block>>, std::vector<std::shared_ptr<block>>> initializeBlocks();
	void sweep();
	//void runFinite();

};
