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
	//cout << singleSiteBasis*basisSize << endl;
	*newSz = identMatrix.kron(Sz1);
	*newSp = identMatrix.kron(Sp1);

	//Transfer the data from the temporary pointers into the permanent variables in the block class
	H=newH;
	Sz=newSz;
	Sp=newSp;

	//We just added a new site, which also increases the basis
	nSites+=1;
	basisSize*=2;

	return;
}

void block::rotateTruncate(matrixReal& transformationMatrix, int maxSize)
{
	int size = min (basisSize, maxSize);
	int singleSiteBasis = 2;
	//cout << singleSiteBasis*size/2 << endl;
	//cout << transformationMatrix.size() << endl;
	//cout << H->size() << endl;
	//cout << transformationMatrix;
	//it's all over 2, b/c the basis size has already been doubled in the enlarge step
	int reducedSize = singleSiteBasis * size/2;
	auto newH = make_shared<matrixReal>(reducedSize,reducedSize);
	*newH = *transformationMatrix.transpose()*(*H*transformationMatrix);
	H = newH; 

	auto newSp = make_shared<matrixReal>(reducedSize,reducedSize);
	*newSp = *transformationMatrix.transpose()*(*Sp*transformationMatrix);
	Sp=newSp;

	auto newSz = make_shared<matrixReal>(reducedSize,reducedSize);
	*newSz = *transformationMatrix.transpose()*(*Sz*transformationMatrix);
	Sz=newSz;

	//*H = *transformationMatrix.transpose()*(*H*transformationMatrix);
	// *Sp = *transformationMatrix.transpose()*(*Sp*transformationMatrix);
	// *Sz = *transformationMatrix.transpose()*(*Sz*transformationMatrix);
	//cout << "Rotated and Truncated H: " << endl << *H;
	//cout << "Rotated and Truncated Sp: " << endl << *Sp;
	//cout << "Rotated and Truncated Sz: " << endl << *Sz;

	basisSize=size;
	
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
	//cout << "System Block Basis Size: "<<sysBlock.basisSize << endl;
	//cout << "Environment Block Basis Size: "<<envBlock.basisSize << endl;

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

shared_ptr<matrixReal> makeReducedDM(matrixReal& groundWfxn)
{	//rho_sys = Tr_env |psi><psi|
	//rho_sys_i,i' = sum_j psi_i,j psi^*_i',j
	//See Schollwock 2011 p.8 or Feguin 2013 p.42

	int nrows = sqrt(groundWfxn.size());
	//cout << groundWfxn.size() << endl;
	auto squarePsi = make_shared<matrixReal>(nrows,nrows);

	for (int ii = 0; ii < nrows; ii++)
	{
		//turn the ground wavefunction into a matrix - it's from the python code
		//it's row major, which may be important? I'm really not sure how or why the python code works
		squarePsi->setSub(0,ii, *groundWfxn.getSub(ii*nrows,0,nrows,1)); //pull out part of the row from groundWfxn
	}
	//cout << endl << "Turn Psi into a matrix: " << endl;
	//cout << *squarePsi;
	//rho = psi * psi^*t
	*squarePsi *= (*squarePsi->transpose());


	return squarePsi; 
}

std::shared_ptr<matrixReal> makeTransformationMatrix(matrixReal& reducedDM, int basisSize, int keepNum)
{
	int	actualKeepNum = min(keepNum, basisSize); //you can't keep more eigenfunctions than exist in the density matrix, only relevant for early steps

	auto transformMatrix = make_shared<matrixReal>(basisSize, actualKeepNum);
	//cout << basisSize << endl << actualKeepNum << endl;
	//now we just need to fill the transformMatrix with elements from the reducedDM, from high eigenvalue to low
	auto extractedColumn = make_shared<matrixReal>(basisSize, 1);
	for (int ii = 0; ii < actualKeepNum; ii++){
		extractedColumn = reducedDM.getSub(0,actualKeepNum-ii-1, basisSize,1); //the high eigenvalues are on the right hand side of the reducedDM
		transformMatrix->setSub(0,ii,*extractedColumn);
	}

	return transformMatrix;
}



double truncationError(std::vector<double>& eigenvals, int basisSize, int keepNum)
{
	double sum = 0;
	int count = 0;
	int offset = 0;
	if (basisSize > keepNum)
	{
		offset=basisSize-keepNum;
	}
	
 	for (auto c : eigenvals){
	    if (count >= offset){
	    	sum+=c;
	    }
	    count ++;
	}

	return 1-sum;
}

int dmrgInfiniteSystem(block& system, int L, int m)
{
	return 0;
}