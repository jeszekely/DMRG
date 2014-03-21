class block{
public:
	int nSites; //chain/number of sites
	int basisSize; //number of m's to keep
	std::shared_ptr<matrixReal> H;
	std::shared_ptr<matrixReal> Sp;
	std::shared_ptr<matrixReal> Sz;

	//Constructor
	block(int, int, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>);
	
	//Copy Constructor
	block(const block&);

};
