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

	auto newH = make_shared<matrixReal>(singleSiteBasis*basisSize, singleSiteBasis*basisSize);

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

shared_ptr<matrixReal> buildSuperblock(block& sysBlock, block& envBlock)
{

	//Step 1, enlarge the system & the environment
	//this is a stupid way of doing things, but it's good enough for now
	//build the initial matrices (H1, Sp1, Sz1) - these are used to expand the hamiltonian, and are not changed
	auto H1 = make_shared<matrixReal>(2,2);

	auto Sp1 = make_shared<matrixReal>(2,2);
	Sp1->element(0,1) = 1.0;

	auto Sz1 = make_shared<matrixReal>(2,2);
	Sz1->element(0,0) = 0.5;
	Sz1->element(1,1) = -0.5;

	//It's quite likely that sysBlock and envBlock will be identical, but let's not force that to be the case

	sysBlock.enlarge(*H1, *Sp1, *Sz1);
	envBlock.enlarge(*H1, *Sp1, *Sz1);
	cout << "System Block Basis Size: "<<sysBlock.basisSize << endl;
	cout << "Environment Block Basis Size: "<<envBlock.basisSize << endl;


	matrixReal sysIdent(sysBlock.basisSize,sysBlock.basisSize);
	sysIdent.makeIdentity();
	
	matrixReal envIdent(envBlock.basisSize,envBlock.basisSize);
	envIdent.makeIdentity();

	//Step 2, make the superblock
	auto superBlock = make_shared<matrixReal>(sysBlock.basisSize*envBlock.basisSize, sysBlock.basisSize*envBlock.basisSize);

	//antiferromagnetic case
	double J = 1;
	double Jz = 1;

	*superBlock = sysBlock.H->kron(envIdent) + envIdent.kron(*envBlock.H) //SysH x I_envSize + I_sysSize x EnvH
					+(sysBlock.Sp->kron(*envBlock.Sp->transpose()))*J/2 //  ( sysSp x envSp*t ) * J/2 
						+(sysBlock.Sp->transpose()->kron(*envBlock.Sp))*J/2 // ( sysSp*t x envSp) * J/2
							+(sysBlock.Sz->kron(*envBlock.Sz))*Jz; // ( sysSz x envSz ) Jz


	//cout << "First 10x10 of Superblock: " << endl << *superBlock << endl;

	return superBlock;
}

int dmrgInfiniteSystem(block& system, int L, int m)
{
	return 0;
}