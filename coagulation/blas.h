#ifndef BLAS_H
#define BLAS_H
#include <omp.h>

#ifdef MKL
#include "mkl.h"
#include "mkl_lapacke.h"
#endif

#ifndef MKL
#include <cblas.h>
#include <lapacke.h>
#endif

#include <complex>
#include <algorithm>

using namespace std;

//Int
void copy(const int & n, const int * x, const int & incx, int * y, const int & incy);

//Double
void orthofactors(const int &, const int &, double * &, double * &, double * &);

void myreducedorthofactors(const double &, const int &, const int &, double * , double * , double * , int &, const int = 0);

void reducedorthofactors(const double &, const int &, const int &, double * &, double * &, double * &, int &);

void my_qr (int m, int n, double * a, int lda, double * matrR, int lda_r );

int gesvd(const int &, const char &, const char &, const lapack_int &, const lapack_int &, double* &, const lapack_int &, double* &, double* &, const lapack_int &, double* &, const lapack_int &, double* &);

void copy(const int &, const double *, const int &, double *, const int &);

void scal(const int &, const double &, double *, const int &);

void axpy(const int &, const double &, const double *, const int &, double *, const int &);

double nrm2(const int &, const double *, const int &);

void gemv(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const int, const int, const double, const double *, const int, const double *, const int, const double, double *, const int);

void gemm(const CBLAS_ORDER, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, const int, const int, const int, const double, const double *, const int, const double *, const int, const double, double *, const int);

double dot_u(const int &, const double *, const int &, const double *, const int &);

double dot_c(const int &, const double *, const int &, const double *, const int &);

double * pseudoinverse(double, int, int, double *);

double imag(const double & );

double real(const double & );

int iamax(const int, const double *, const int);

int gesv( int , int , int nrhs, double* a, int lda, int* ipiv, double* b, int ldb );

void ger (const CBLAS_ORDER, const int, const int , const double, const double *, const int , const double *, const int, double *, const int );

void transpose(const double *, const int, const int, double *, int = 0, int = 0);

#endif
