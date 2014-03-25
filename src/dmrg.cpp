#include <iostream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

#include "matrix.hpp"
#include "block.hpp"
#include "dmrg.hpp"
#include "utilities.hpp"

using namespace std;

// vector <shared_ptr<block>> Blocks(nBlocks);
// for (auto &blk : Blocks)
// 	blk = make_shared<block>(block_args);


infiniteSystem::infiniteSystem(int sl, vector<int> kn) : sysLength(sl) , sweepSizes(kn) {}


//***************************************************************************
//**********************	Infinite System 	*****************************
//***************************************************************************


void infiniteSystem::runInfinite()
{
	double energy, error;
	int maxEigenStates=sweepSizes.front();

	  shared_ptr<matrixReal> H, Sp, Sz, eH, eSp, eSz;
	  tie(H, Sp, Sz) = makeSingleSiteOperators();
	  //Might need these separate, not sure yet
	  tie(eH, eSp, eSz) = makeSingleSiteOperators();

	  //Create the blocks
	  block sysBlock(1, 2, H, Sp, Sz);
	  block envBlock(1, 2, eH, eSp, eSz);

	for (int currentLen = 0; 2*(currentLen+1) < sysLength; currentLen++){
		cout << "System Length = " << 2*sysBlock.nSites+2 << endl;

		tie(error,energy)=infiniteDMRGStep(sysBlock, envBlock, maxEigenStates);

		cout << "Truncation Error : " << error << endl;
		cout << "E/L = " << std::setprecision(10) << energy/(sysBlock.nSites*2) << endl;
		//sysBlock.H->printMem();
		cout << endl;
	}
  return;
}

std::tuple<double, double> infiniteSystem::infiniteDMRGStep(block& sysBlock, block& envBlock, int maxEigenStates)
{	//Propagate the system for one step

	//*** Build the Superblock and get the ground state energy and wavefunction***
	double energy, error;



 	enlarge(sysBlock);
    enlarge(envBlock);

    shared_ptr<matrixReal> groundState;
    tie(energy,groundState) = buildSuperblock(sysBlock, envBlock);

    shared_ptr<matrixReal> reducedDM;
    tie(error,reducedDM)= makeReducedDM(sysBlock.basisSize, envBlock.basisSize, *groundState,maxEigenStates);

    //***Diagonalize Reduced Density Matrix***
    //***Make transformation matrix from reduced density matrix by keeping only up to maxKeepNum eigenstates.***
    makeTransformationMatrix(sysBlock, *reducedDM,maxEigenStates);
    envBlock = sysBlock; //Since we don't rotate/truncate the environment block, we just overwrite it with the sysblock


	return make_tuple(error, energy);
}

std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>> infiniteSystem::makeSingleSiteOperators()
{
	auto H1 = make_shared<matrixReal>(2,2);

	auto Sp1 = make_shared<matrixReal>(2,2);
	Sp1->element(0,1) = 1.0;

	auto Sz1 = make_shared<matrixReal>(2,2);
	Sz1->element(0,0) = 0.5;
	Sz1->element(1,1) = -0.5;

	return make_tuple(H1, Sp1, Sz1);
}


void infiniteSystem::enlarge(block& b)
{
    shared_ptr<matrixReal> H1, Sp1, Sz1;
    tie(H1, Sp1, Sz1) = makeSingleSiteOperators();


	//H1, Sz1, and Sp1 all refer to the original, 2x2 matrices. They should never change.
	//construct identity matrix for upscaling stuff
	matrixReal identMatrix(b.basisSize,b.basisSize);
	identMatrix.makeIdentity();

	//construct 2x2 identity matrix - 2x2 is only appropriate for the spin system, change later for other systems
	int singleSiteBasis = 2;
	matrixReal smallIdent(singleSiteBasis,singleSiteBasis);
	smallIdent.makeIdentity();


	//antiferromagnetic case
	double J = 1;
	double Jz = 1;

	auto newH = make_shared<matrixReal>(singleSiteBasis*b.basisSize, singleSiteBasis*b.basisSize);
	//See Chapter 2 The Density Matrix Renormalization Group (Adrian E. Feiguin) eq. 2.6, 2.8, 2.11
	*newH = b.H->kron(smallIdent) //H x I_2
						+ identMatrix.kron(*H1)  //I_basisSize x H_1
							+ (b.Sp->kron(*(Sp1->transpose())) //Sp x Sp1^(*t)
								+ (b.Sp->transpose())->kron(*Sp1))*(J/2.) //Sp^(*t) x Sp1
									+(b.Sz->kron(*Sz1))*Jz; //Sz x Sz1 * Jz

	//Scale up Sz and Sp to the proper size
	auto newSz = make_shared<matrixReal>(singleSiteBasis*b.basisSize, singleSiteBasis*b.basisSize);
	auto newSp = make_shared<matrixReal>(singleSiteBasis*b.basisSize, singleSiteBasis*b.basisSize);

	*newSz = identMatrix.kron(*Sz1);
	*newSp = identMatrix.kron(*Sp1);

	//Transfer the data from the temporary pointers into the permanent variables in the block class
	b.H=newH;
	b.Sz=newSz;
	b.Sp=newSp;

	//We just added a new site, which also increases the basis
	b.nSites+=1;
	b.basisSize*=2;

	return;
}


std::tuple<double, std::shared_ptr<matrixReal>> infiniteSystem::buildSuperblock(block& enlargedSysBlock, block& enlargedEnvBlock)
{

	matrixReal sysIdent(enlargedSysBlock.basisSize,enlargedSysBlock.basisSize);
	sysIdent.makeIdentity();

	matrixReal envIdent(enlargedEnvBlock.basisSize,enlargedEnvBlock.basisSize);
	envIdent.makeIdentity();

	auto superBlock = make_shared<matrixReal>(enlargedSysBlock.basisSize*enlargedEnvBlock.basisSize, enlargedSysBlock.basisSize*enlargedEnvBlock.basisSize);

	//antiferromagnetic case
	double J = 1;
	double Jz = 1;

	*superBlock = enlargedSysBlock.H->kron(envIdent) + sysIdent.kron(*enlargedEnvBlock.H) //SysH x I_envSize + I_sysSize x EnvH
					+(enlargedSysBlock.Sp->kron(*enlargedEnvBlock.Sp->transpose()))*J/2 //  ( sysSp x envSp*t ) * J/2
						+(enlargedSysBlock.Sp->transpose()->kron(*enlargedEnvBlock.Sp))*J/2 // ( sysSp*t x envSp) * J/2
							+(enlargedSysBlock.Sz->kron(*enlargedEnvBlock.Sz))*Jz; // ( sysSz x envSz ) Jz

	// *** Diagonalize the Superblock and get the eigenvalues ***
	//vector <double> superBlockVals(enlargedSysBlock.basisSize*enlargedEnvBlock.basisSize,0.0);
	vector <double> superBlockVals(1,0.0);
	//superBlock->diagonalize(superBlockVals.data());

	superBlock->diagonalize(superBlockVals.data(),1,1);

	// *** Get the ground state of the diagonalized superblock
	auto groundState = make_shared<matrixReal>(*superBlock->getSub(0,0,enlargedSysBlock.basisSize*enlargedEnvBlock.basisSize,1));
	// ** Return the ground state energy and the ground state itself ***
	return make_tuple(*superBlockVals.data(), groundState);
}



std::tuple<double, std::shared_ptr<matrixReal>> infiniteSystem::makeReducedDM(int sysBasisSize, int envBasisSize, matrixReal& groundWfxn, int maxEigenStates)
{
	// *** Make Reduced Density Matrix ***

    const int environment_basis = groundWfxn.size()/sysBasisSize;
    assert (environment_basis == envBasisSize);


	auto squarePsi = make_shared<matrixReal>(sysBasisSize,sysBasisSize);
    //We transpose the first matrix. Basically, the 2nd index is the one that moves first (a consequence of the way the kronecker product is written)
	dgemm_("T", "N", sysBasisSize, sysBasisSize, envBasisSize, 1.0, groundWfxn.data(), envBasisSize, groundWfxn.data(), envBasisSize, 0.0, squarePsi->data(), sysBasisSize);

	// *** Diagonalize the Reduced DM, and get the eigenvalues ***
    vector <double> reducedDMVals(sysBasisSize,0.0);
    squarePsi->diagonalize(reducedDMVals.data());

    // *** Calculate the truncation error ***
    double error = accumulate(reducedDMVals.begin(),reducedDMVals.begin()+max(sysBasisSize-maxEigenStates,0),0.0);

	return make_tuple(error, squarePsi);
}



void infiniteSystem::makeTransformationMatrix(block& sysBlock, matrixReal& reducedDM, int maxEigenStates)
{
	// *** Make the Transformation Matrix ***
	int	actualKeepNum = min(maxEigenStates, sysBlock.basisSize);

	// *** Extract the transformation matrix from the right hand side of the reduced DM
	auto transformationMatrix = reducedDM.getSub(0,sysBlock.basisSize-actualKeepNum,sysBlock.basisSize,actualKeepNum);

	sysBlock.H = make_shared<matrixReal>(*transformationMatrix | (*sysBlock.H * *transformationMatrix));
	sysBlock.Sp = make_shared<matrixReal>(*transformationMatrix | (*sysBlock.Sp * *transformationMatrix));
	sysBlock.Sz = make_shared<matrixReal>(*transformationMatrix | (*sysBlock.Sz * *transformationMatrix));
	sysBlock.basisSize=actualKeepNum;


	return;
}

void infiniteSystem::graphics(int sysBlockLength, int envBlockLength, bool direction)
{
	string DMRGState;
	DMRGState.append(size_t(sysBlockLength),'=');
	DMRGState.append("**");
	DMRGState.append(size_t(envBlockLength),'-');
	if (!direction) DMRGState = string(DMRGState.rbegin(),DMRGState.rend());
	cout << DMRGState << endl;
}

//***************************************************************************
//**********************	Finite System 	*********************************
//***************************************************************************

finiteSystem::finiteSystem(int sl, vector<int> kn) : infiniteSystem(sl, kn) {}


std::tuple< vector<std::shared_ptr<block>>, vector<std::shared_ptr<block>>> finiteSystem::initializeBlocks()
{

	shared_ptr<matrixReal> H, Sp, Sz, eH, eSp, eSz;
 	tie(H, Sp, Sz) = makeSingleSiteOperators();
	//Might need these separate, not sure yet
	tie(eH, eSp, eSz) = makeSingleSiteOperators();

	double energy, error;

	vector <shared_ptr<block>> leftBlocks(sysLength);
 	vector <shared_ptr<block>> rightBlocks(sysLength);

 	leftBlocks[0] = make_shared<block>(1, 2, H, Sp, Sz);
 	rightBlocks[0] = make_shared<block>(1, 2, eH, eSp, eSz);

 	cout << "Creating Initial " << sysLength << " blocks." << endl;
    for (int currentLen = 0; 2*(currentLen+1) < sysLength; currentLen++){
    	// ** Copy the current blocks into new variables
    	block newLeftBlock(*leftBlocks[currentLen]);
	 	block newRightBlock(*rightBlocks[currentLen]);

	 	cout << "System Length = " << 2*newLeftBlock.nSites+2 << endl;
	    graphics(newLeftBlock.nSites, newRightBlock.nSites, 1);
	    //** This expands both the left and right block, which we need to assign to
	    tie(error,energy)=infiniteDMRGStep(newLeftBlock, newRightBlock, sweepSizes.front()); //the initial step currently uses the first value in sweepSizes

	    cout << "Truncation Error : " << error << endl;
	    cout << "E/L = " << std::setprecision(10) << energy/(newLeftBlock.nSites*2) << endl;
	    //newLeftBlock.H->printMem();
	    cout << endl;
	    leftBlocks[currentLen+1] = make_shared<block>(newLeftBlock);
	    rightBlocks[currentLen+1] = make_shared<block>(newLeftBlock); //at this stage, the two sides are identical
    }
    cout << "Initial Chain Complete" << endl << endl << endl;

	vector<shared_ptr<block>> *sysList = &leftBlocks;
	vector<shared_ptr<block>> *envList = &rightBlocks;

	return make_tuple(*sysList, *envList);



}
void finiteSystem::sweep()
{
	// ** Pull out and copy the biggest system block element
	//Set the intitial system block

	vector<shared_ptr<block>> sysList;
	vector<shared_ptr<block>> envList;

	tie(sysList, envList)=initializeBlocks();

	block sysBlock(*sysList.at(sysLength/2-1));

	double energy, error;

	for (auto maxEigenStates : sweepSizes){
		cout << "Beginning Sweep With Max Basis Size: " << maxEigenStates << endl << endl;
		bool direction = 1;
	    for (int cycle = 0; cycle<(2*sysLength-8); cycle++) //One complete sweep
	    {
	    	//get the envBlock corresponding to the current sysBlock
	    	block envBlock(*envList.at(sysLength-sysBlock.nSites-3)); // -2 for the 2 active sites, -1 for the 0 based array we use


	    	if (envBlock.nSites == 1)
	    	{ //swap
	    		std::swap(sysList, envList);
	    		std::swap(sysBlock, envBlock);
	    		direction = !direction;
	    	}


			cout << "System Block Size: "<< sysBlock.nSites << ", ";
	    	cout << "Environment Block Size: "<< envBlock.nSites << endl;
   		    graphics(sysBlock.nSites, envBlock.nSites, direction);

	    	tie(error,energy)=infiniteDMRGStep(sysBlock, envBlock, maxEigenStates);

			cout << "Truncation Error : " << error << endl;
		    cout << "E/L = " << std::setprecision(10) << energy/sysLength << endl;
		    sysBlock.H->printMem();
	    	sysList.at(sysBlock.nSites-1) = make_shared<block>(sysBlock);
		    cout << endl;

	    }
	}
}

