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

//Constructor
block::block(int ns, int bs, std::shared_ptr<matrixReal> h, std::shared_ptr<matrixReal> sp, std::shared_ptr<matrixReal> sz) : nSites(ns), basisSize(bs), H(h), Sp(sp), Sz(sz) {}

//Copy Constructor
block::block(const block& b) : nSites(b.nSites), basisSize(b.basisSize), H(b.H), Sp(b.Sp), Sz(b.Sz)  {}
