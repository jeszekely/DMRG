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
#include "utilities.hpp"

using namespace std;

block::block(int ns, int bs, int kn, std::shared_ptr<matrixReal> h, std::shared_ptr<matrixReal> sp, std::shared_ptr<matrixReal> sz) : nSites(ns), basisSize(bs), maxKeepNum(kn), H(h), Sp(sp), Sz(sz) {}

std::tuple<double, double> block::infiniteDMRGStep(block& environmentBlock)
{	//Propagate the system for one step

	//*** Build the Superblock and get the ground state energy and wavefunction***
	shared_ptr<matrixReal> groundState;
	double energy;
	tie(energy,groundState) = this->buildSuperblock(environmentBlock);
	


    //*** Make Reduced Density Matrix, then Diagonalize it, and get the truncation error while we're at it***
   	shared_ptr<matrixReal> reducedDM;
	double error;
   	tie(error,reducedDM)= this->makeReducedDM(*groundState);

    //***Diagonalize Reduced Density Matrix***
 
    //***Make transformation matrix from reduced density matrix by keeping only up to maxKeepNum eigenstates.***
    this->makeTransformationMatrix(*reducedDM, environmentBlock);

    
	return make_tuple(error, energy); 
}

std::tuple<double, std::shared_ptr<matrixReal>> block::buildSuperblock(block& envBlock)
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
	this->enlarge(*H1, *Sp1, *Sz1);
	envBlock.enlarge(*H1, *Sp1, *Sz1);


	matrixReal sysIdent(basisSize,basisSize);
	sysIdent.makeIdentity();

	matrixReal envIdent(envBlock.basisSize,envBlock.basisSize);
	envIdent.makeIdentity();

	//Step 2, make the superblock
	auto superBlock = make_shared<matrixReal>(basisSize*envBlock.basisSize, basisSize*envBlock.basisSize);

	//antiferromagnetic case
	double J = 1;
	double Jz = 1;
	*superBlock = H->kron(envIdent) + envIdent.kron(*envBlock.H) //SysH x I_envSize + I_sysSize x EnvH
					+(Sp->kron(*envBlock.Sp->transpose()))*J/2 //  ( sysSp x envSp*t ) * J/2
						+(Sp->transpose()->kron(*envBlock.Sp))*J/2 // ( sysSp*t x envSp) * J/2
							+(Sz->kron(*envBlock.Sz))*Jz; // ( sysSz x envSz ) Jz


	
	// *** Diagonalize the Superblock and get the eigenvalues ***
	vector <double> superBlockVals(basisSize*basisSize,0.0);
	superBlock->diagonalize(superBlockVals.data());

	// *** Get the ground state of the diagonalized superblock
	auto groundState = make_shared<matrixReal>(*superBlock->getSub(0,0,basisSize*basisSize,1));

	// ** Return the ground state energy and the ground state itself ***
	return make_tuple(*superBlockVals.data(), groundState);
}

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


std::tuple<double, std::shared_ptr<matrixReal>> block::makeReducedDM(matrixReal& groundWfxn)
{
	// *** Make Reduced Density Matrix ***
	const int nrows = sqrt(groundWfxn.size());
	auto squarePsi = make_shared<matrixReal>(nrows,nrows);
	dgemm_("N", "T", nrows, nrows, nrows, 1.0, groundWfxn.data(), nrows, groundWfxn.data(), nrows, 0.0, squarePsi->data(), nrows);
	
	// *** Diagonalize the Reduced DM, and get the eigenvalues ***
    vector <double> reducedDMVals(basisSize,0.0);
    squarePsi->diagonalize(reducedDMVals.data());

    // *** Calculate the truncation error ***
    double error = accumulate(reducedDMVals.begin(),reducedDMVals.begin()+max(basisSize-maxKeepNum,0),0.0);

	return make_tuple(error, squarePsi);
}


void block::makeTransformationMatrix(matrixReal& reducedDM, block& environmentBlock)
{	
	// *** Make the Transformation Matrix ***
	int	actualKeepNum = min(maxKeepNum, basisSize);
	auto transformationMatrix = reducedDM.getSub(0,basisSize-actualKeepNum,basisSize,actualKeepNum);
	
	// *** While we're here, rotate and truncate the system and environment blocks ***
	this->rotateTruncate(*transformationMatrix);
    environmentBlock.rotateTruncate(*transformationMatrix);

	return; 
}


void block::rotateTruncate(matrixReal& transformationMatrix)
{
	int size = min (basisSize, maxKeepNum);
	H = make_shared<matrixReal>(transformationMatrix | *H * transformationMatrix);
	Sp = make_shared<matrixReal>(transformationMatrix | *Sp * transformationMatrix);
	Sz = make_shared<matrixReal>(transformationMatrix | *Sz * transformationMatrix);
	basisSize=size;
	return;
}

<<<<<<< HEAD
=======
double truncationError(std::vector<double>& eigenvals, int basisSize, int keepNum)
{
	//This probably doesn't need to be a separate function
	return accumulate(eigenvals.begin(),eigenvals.begin()+max(basisSize-keepNum,0),0.0);
}
>>>>>>> 59b5cd9e00105a6990c0d0dae2e528a15b9691a7

