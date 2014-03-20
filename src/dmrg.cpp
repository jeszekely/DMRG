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

// vector <shared_ptr<block>> Blocks(nBlocks);
// for (auto &blk : Blocks)
// 	blk = make_shared<block>(block_args);


block::block(int ns, int bs, int kn, std::shared_ptr<matrixReal> h, std::shared_ptr<matrixReal> sp, std::shared_ptr<matrixReal> sz) : nSites(ns), basisSize(bs), maxKeepNum(kn), H(h), Sp(sp), Sz(sz) {}
block::block(const block& b) : nSites(b.nSites), basisSize(b.basisSize), maxKeepNum(b.maxKeepNum), H(b.H), Sp(b.Sp), Sz(b.Sz)  {}

finiteSystem::finiteSystem(int L, int kn, vector<std::shared_ptr<block>> lb, vector<std::shared_ptr<block>> rb) : sysLength(L), maxKeepNum(kn), leftBlocks(lb), rightBlocks(rb) {}


void finiteSystem::initializeBlocks()
{
	// ** Create Initial Blocks **
	std::shared_ptr<matrixReal> H1(new matrixReal(2,2));
	auto rH1 = make_shared<matrixReal>(*H1);

	std::shared_ptr<matrixReal> Sp1(new matrixReal(2,2));
	Sp1->element(0,1) = 1.0;
	auto rSp1 = make_shared<matrixReal>(*Sp1);

	std::shared_ptr<matrixReal> Sz1(new matrixReal(2,2));
	Sz1->element(0,0) = 0.5;
	Sz1->element(1,1) = -0.5;
	auto rSz1 = make_shared<matrixReal>(*Sz1);

	// ** Create the initial blocks of length 1 **
 	double energy, error;
 	leftBlocks[0] = make_shared<block>(1, 2, maxKeepNum, H1, Sp1, Sz1);
 	rightBlocks[0] = make_shared<block>(1, 2, maxKeepNum, rH1, rSp1, rSz1);

 	// ** Create the initital, parallel system **
 	cout << "Creating Initial " << sysLength << " blocks." << endl;
    for (int currentLen = 0; 2*(currentLen+1) < sysLength; currentLen++){
    	// ** Copy the current blocks into new variables
    	block newLeftBlock(*leftBlocks[currentLen]);
	 	block newRightBlock(*rightBlocks[currentLen]);

	 	cout << "System Length = " << 2*newLeftBlock.nSites+2 << endl;
	    
	    //** This expands both the left and right block, which we need to assign to 
	    tie(error,energy)=newLeftBlock.infiniteDMRGStep(newRightBlock);
	    
	    cout << "Truncation Error : " << error << endl;
	    cout << "E/L = " << std::setprecision(10) << energy/(newLeftBlock.nSites*2) << endl;
	    newLeftBlock.H->printMem();
	    cout << "******************************************************************************" << endl;
	    leftBlocks[currentLen+1] = make_shared<block>(newLeftBlock);
	    rightBlocks[currentLen+1] = make_shared<block>(newLeftBlock); //at this stage, the two sides are identical
    }


    cout << "Initial Chain Complete" << endl << endl << endl;


	vector<shared_ptr<block>> *sysList = &leftBlocks;
	vector<shared_ptr<block>> *envList = &rightBlocks;

    // ** Pull out and copy the biggest system block element
	//Set the intitial system block
	block sysBlock(*sysList->at(sysLength/2-1));



    for (int cycle = 0; cycle<(2*sysLength-8); cycle++) //One complete sweep
    {
    	//get the envBlock corresponding to the current sysBlock
    	block envBlock(*envList->at(sysLength-sysBlock.nSites-3)); // -2 for the 2 active sites, -1 for the 0 based array we use
    	
   		
    	if (envBlock.nSites == 1)
    	{ //swap
    		std::swap(sysList, envList);
    		std::swap(sysBlock, envBlock);
    	}

		cout << "System Block Size: "<< sysBlock.nSites << endl;
    	cout << "Environment Block Size: "<< envBlock.nSites << endl;

    	tie(error,energy)=sysBlock.infiniteDMRGStep(envBlock);

		cout << "Truncation Error : " << error << endl;
	    cout << "E/L = " << std::setprecision(10) << energy/sysLength << endl;
	    //newLeftBlock.H->printMem();
    	sysList->at(sysBlock.nSites-1) = make_shared<block>(sysBlock);
    	//cout << *(sysBlock.H);
	    cout << "******************************************************************************" << endl;

	   	//cout << cycle << endl;
    }



	// for (auto &blk : leftBlocks){
 	// 		blk = make_shared<block>(1, 2, maxKeepNum, H1, Sp1, Sz1);
 	// 	}

	return;
}

std::tuple<double, double> block::infiniteDMRGStep(block& environmentBlock)
{	//Propagate the system for one step

	//*** Build the Superblock and get the ground state energy and wavefunction***

	shared_ptr<matrixReal> groundState;
	double energy;
	tie(energy,groundState) = this->buildSuperblock(environmentBlock);
	//cout << "Ground State: " << endl << *groundState;
	//*** Make Reduced Density Matrix, then Diagonalize it, and get the truncation error while we're at it***
   	shared_ptr<matrixReal> reducedDM;
	double error;
   	tie(error,reducedDM)= this->makeReducedDM(*groundState);
   	//cout << *reducedDM;
    //***Diagonalize Reduced Density Matrix***
    //***Make transformation matrix from reduced density matrix by keeping only up to maxKeepNum eigenstates.***
    this->makeTransformationMatrix(*reducedDM);

    
	return make_tuple(error, energy); 
}

std::tuple<double, std::shared_ptr<matrixReal>> block::buildSuperblock(block& envBlock)
{
	cout << "Build Superblock" << endl;
	//Step 1, enlarge the system & the environment
	//this is a stupid way of doing things, but it's good enough for now
	//build the initial matrices (H1, Sp1, Sz1) - these are used to expand the hamiltonian, and are not changed
	auto H1 = make_shared<matrixReal>(2,2);

	auto Sp1 = make_shared<matrixReal>(2,2);
	Sp1->element(0,1) = 1.0;

	auto Sz1 = make_shared<matrixReal>(2,2);
	Sz1->element(0,0) = 0.5;
	Sz1->element(1,1) = -0.5;

	// cout << "SysBlock (Not Enlarged)" << endl;
	// cout << H->size() << endl;
	// cout << *H << endl;
	this->enlarge(*H1, *Sp1, *Sz1);
	// cout << *H << endl;
	// cout << *envBlock.H << endl;

	envBlock.enlarge(*H1, *Sp1, *Sz1);
	
	// cout << "SysBlock (Enlarged)" << endl;
	// cout << H->size() << endl;
	// cout << *H << endl;

	// cout << "EnvBlock (Enlarged)" << endl;
	// cout << envBlock.H->size() << endl;
	// cout << *envBlock.H << endl;


	matrixReal sysIdent(basisSize,basisSize);
	sysIdent.makeIdentity();

	matrixReal envIdent(envBlock.basisSize,envBlock.basisSize);
	envIdent.makeIdentity();

	//Step 2, make the superblock
	auto superBlock = make_shared<matrixReal>(basisSize*envBlock.basisSize, basisSize*envBlock.basisSize);
	//antiferromagnetic case
	double J = 1;
	double Jz = 1;



	*superBlock = H->kron(envIdent) + sysIdent.kron(*envBlock.H) //SysH x I_envSize + I_sysSize x EnvH
					+(Sp->kron(*envBlock.Sp->transpose()))*J/2 //  ( sysSp x envSp*t ) * J/2
						+(Sp->transpose()->kron(*envBlock.Sp))*J/2 // ( sysSp*t x envSp) * J/2
							+(Sz->kron(*envBlock.Sz))*Jz; // ( sysSz x envSz ) Jz

	// *** Diagonalize the Superblock and get the eigenvalues ***
	vector <double> superBlockVals(basisSize*envBlock.basisSize,0.0);
    //cout << "Superblock" << endl;
    //cout << *superBlock;
	superBlock->diagonalize(superBlockVals.data());
    //cout << "End Diag" << endl;

	// *** Get the ground state of the diagonalized superblock
	//cout <<"GroundWfxn  " << basisSize << "\t" << envBlock.basisSize << "\t" << basisSize*envBlock.basisSize << endl;
	auto groundState = make_shared<matrixReal>(*superBlock->getSub(0,0,basisSize*envBlock.basisSize,1));
	// ** Return the ground state energy and the ground state itself ***
	return make_tuple(*superBlockVals.data(), groundState);
}

void block::enlarge(matrixReal &H1, matrixReal &Sp1, matrixReal &Sz1)
{
	//H1, Sz1, and Sp1 all refer to the original, 2x2 matrices. They should never change.
	cout << "Enlarge Block" << endl;
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
	// cout << *H;
	// cout << (H->kron(smallIdent));
	// cout << (identMatrix.kron(H1));
	// cout << (Sp->kron(*Sp1.transpose()));
	// cout << ((Sp->transpose())->kron(Sp1));
	// cout << (Sz->kron(Sz1));

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
    cout << "Make Reduced DM" <<endl;

	//cout << basisSize << "," << groundWfxn.size() << "," <<  groundWfxn.size()/basisSize << endl;
	auto squarePsi = make_shared<matrixReal>(basisSize,basisSize);
	matrixReal ground_copy(groundWfxn);
	dgemm_("N", "T", basisSize, basisSize, basisSize, 1.0, groundWfxn.data(), basisSize, ground_copy.data(), basisSize, 0.0, squarePsi->data(), basisSize);
	
	//cout << "Square Psi: " << endl << *squarePsi;
	// *** Diagonalize the Reduced DM, and get the eigenvalues ***
    vector <double> reducedDMVals(basisSize,0.0);
    squarePsi->diagonalize(reducedDMVals.data());


    //cout << *squarePsi;
 //    for (auto c : reducedDMVals)
 //    std::cout << c << ' ';
	// cout << endl;


    // *** Calculate the truncation error ***
    double error = accumulate(reducedDMVals.begin(),reducedDMVals.begin()+max(basisSize-maxKeepNum,0),0.0);

	return make_tuple(error, squarePsi);
}


void block::makeTransformationMatrix(matrixReal& reducedDM)
{	
	// *** Make the Transformation Matrix ***
	int	actualKeepNum = min(maxKeepNum, basisSize);
	
	// *** Extract the transformation matrix from the right hand side of the reduced DM
	auto transformationMatrix = reducedDM.getSub(0,basisSize-actualKeepNum,basisSize,actualKeepNum);

	// *** While we're here, rotate and truncate the system block ***
	this->rotateTruncate(*transformationMatrix);


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


