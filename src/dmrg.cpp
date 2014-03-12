#include <iostream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"
#include "dmrg.hpp"

using namespace std;

block::block(int ns, int bs, std::shared_ptr<matrixReal> h, std::shared_ptr<matrixReal> sp, std::shared_ptr<matrixReal> sz) : nSites(ns), basisSize(bs), H(h), Sp(sp), Sz(sz) {}

void block::enlarge(matrixReal &H1, matrixReal &Sp1, matrixReal &Sz1)
{
	//H1, Sz1, and Sp1 all refer to the original, 2x2 matrices. They should never change.

	//construct identity matrix for upscaling stuff
	matrixReal identMatrix(basisSize,basisSize);
	identMatrix.makeIdentity();

	//construct 2x2 identity matrix - 2x2 is only appropriate for the spin system, change later for other systems
	int singleSiteBasis = 2;
	matrixReal smallIdent(singleSiteBasis,singleSiteBasis);
	smallIdent.makeIdentity();


	//antiferromagnetic case
	double J = 1;
	double Jz = 1;

	std::shared_ptr<matrixReal> newH(new matrixReal (singleSiteBasis*basisSize, singleSiteBasis*basisSize));
	//See Chapter 2 The Density Matrix Renormalization Group (Adrian E. Feiguin) eq. 2.6, 2.8, 2.11
	*newH = H->kron(smallIdent) //H x I_2
						+ identMatrix.kron(H1)  //I_basisSize x H_1
							+ (Sp->kron(*Sp1.transpose()) //Sp x Sp1^(*t)
								+ (Sp->transpose())->kron(Sp1))*(J/2.) //Sp^(*t) x Sp1
									+(Sz->kron(Sz1))*Jz; //Sz x Sz1 * Jz

	//Scale up Sz and Sp to the proper size
	auto newSz = make_shared<matrixReal>(singleSiteBasis*basisSize, singleSiteBasis*basisSize);
	auto newSp = make_shared<matrixReal>(singleSiteBasis*basisSize, singleSiteBasis*basisSize);

	*newSz = identMatrix.kron(*Sz);
	*newSp = identMatrix.kron(*Sp);

	//Transfer the data from the temporary pointers into the permanent variables in the block class
	H=newH;
	Sz=newSz;
	Sp=newSp;

	//We just added a new site, which also increases the basis
	nSites+=1;
	basisSize*=2;

	return;
}

int dmrgInfiniteSystem(block& system, int L, int m)
{
	return 0;
}