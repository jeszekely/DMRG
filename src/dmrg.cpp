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

block::block(int ns, int bs, std::shared_ptr<matrixReal> h, std::shared_ptr<matrixReal> ju, std::shared_ptr<matrixReal> jl) : nSites(ns), basisSize(bs), H(h), JU(ju), JL(jl) {}

void block::enlarge(block &o)
{
	return;
}

int dmrgInfiniteSystem(block& system, int L, int m)
{
	return 0;
}