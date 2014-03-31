#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "davidson.hpp"

using namespace std;

genMatrix::genMatrix(size_t nr, size_t nc, std::function<vectorMatrix(vectorMatrix&)> h, vector<double> &Diag) : nrows(nr), ncols(nc), H(h),diags(Diag) {}

double genMatrix::diagElem(int ii) {return diags[ii];}

int genMatrix::nc() {return ncols;}

int genMatrix::nr() {return nrows;}

vectorMatrix genMatrix::operator*(vectorMatrix& o){return H(o);}

Davidson::Davidson(genMatrix h, int iVecs, int nVecs, int maxIts, double err) : H(h), initialVecs(iVecs), numVecs(nVecs), maxIterations(maxIts), tolerance(err) {}

tuple<std::shared_ptr<vectorMatrix>,vector<double>> Davidson::diagonalize()
{
  const int n = H.nr(); //Size of a vector

  //eigenvector and eigenvalue return arrays
  auto eigVecs = make_shared<vectorMatrix>(n,numVecs);
  vector<double> eigVals(numVecs, 0.0);

  //initial guess set to unit vectors
  auto trials = make_shared<vectorMatrix>(n, initialVecs);
  trials->makeIdentity();

  for (int iter = 0; iter < maxIterations; ++iter) {

    //Apply matrix to guess and make subspace matrix, diagonalize
    vectorMatrix Ht(H * *trials);
    vectorMatrix Hsub( *trials | Ht );
    const int subsize = Hsub.nc();
    vector<double> eigs(subsize, 0.0);
    Hsub.diagonalize(eigs.data());

    vectorMatrix psi( *trials*Hsub );    // current best guesses for eigenvectors
    vectorMatrix sigma( Ht*Hsub ); // transformed sigma vectors

    copy_n(psi.data(), n*numVecs, eigVecs->data());
    copy_n(eigs.data(), numVecs, eigVals.data());

    vector<shared_ptr<matrixReal>> new_trial_vectors;

    for (int ii = 0; ii < numVecs; ++ii) {
      //matrixReal residual(n, 1);
      //daxpy_(n, -eigs[ii], &psi(0,ii), 1, &residual(0,0), 1);
      //daxpy_(n, 1.0, &sigma(0,ii), 1, &residual(0,0), 1);
      matrixReal residual( (psi.vec(ii)*eigs[ii]) - sigma.vec(ii) );
      cout << residual;

      const double residual_norm = residual.variance();

      cout << setw(6) << iter << setw(6) << ii
                              << fixed << setw(22) << setprecision(12) << eigs[ii]
                              << scientific << setw(16) << setprecision(8) << residual_norm << endl;

      if ( residual_norm > tolerance) {
        auto trial_vector = make_shared<matrixReal>(n, 1);

        // form new guess
        const double en = eigs[ii];
        for (int kk = 0; kk < n; ++kk)
          trial_vector->element(kk, 0) = residual(kk,0) / min(en - H.diagElem(kk), -0.1);

        new_trial_vectors.push_back(trial_vector);
      }
    }
    if (numVecs != 0) cout << endl;

    if (new_trial_vectors.empty()) {
      cout << "Converged!" << endl;
      break;
    }

    auto new_trials = make_shared<vectorMatrix>(n, trials->nc() + new_trial_vectors.size());
    copy_n(trials->data(), trials->size(), new_trials->data());
    for (int ii = 0; ii < new_trial_vectors.size(); ++ii)
      copy_n(new_trial_vectors[ii]->data(), new_trial_vectors[ii]->nr(), &new_trials->element(0, trials->nc() + ii));

    trials = new_trials->canonical_orthogonalization();
  }
  return make_tuple(eigVecs,eigVals);
}
