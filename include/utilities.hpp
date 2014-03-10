#include <complex>
#include <memory>

//BLAS
extern "C"
{
	void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
		const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
		const double* beta, double* c, const int* ldc);

  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
}

//DMRG interface, shamelessly lifted from BAGEL's src/util/f77.h file
namespace
{
	void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
             const double alpha, const double* a, const int lda, const double* b, const int ldb,
             const double beta, double* c, const int ldc) { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }
	void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
             const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb,
             const double beta, std::unique_ptr<double []>& c, const int ldc)
             { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

	void dsyev_(const char* a, const char* b, const int c, double* d, const int e, double* f, double* g, const int h, int& i)
             { ::dsyev_(a,b,&c,d,&e,f,g,&h,&i);}
	void dsyev_(const char* a, const char* b, const int c, std::unique_ptr<double []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<double []>& g, const int h, int& i)
             { ::dsyev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,&i);}
}