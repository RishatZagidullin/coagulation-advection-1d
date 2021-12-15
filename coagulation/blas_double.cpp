#include "blas.h"
#include <iostream>
using namespace std;

void orthofactors(const int & m, const int &n, double * &a, double * &U, double * &VT)
{
    double *sigma, *superb;
    if (m <= n)
    {
        U = (double *) malloc(m * m * sizeof(double));
        VT = (double *) malloc(m * n * sizeof(double));
        sigma = (double *) malloc(m * sizeof(double));
        superb = (double *) malloc((m - 1) * sizeof(double));
        gesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, a, m, sigma, U, m, VT, m, superb);
        for (int i = 0; i < m; i++)
        {
            scal(m, sigma[i], U + i * m, 1);
        }
    }
    else
    {
        U = (double *) malloc(m * n * sizeof(double));
        VT = (double *) malloc(n * n * sizeof(double));
        sigma = (double *) malloc(n * sizeof(double));
        superb = (double *) malloc((n - 1) * sizeof(double));
        gesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, a, m, sigma, U, m, VT, m, superb);
        for (int i = 0; i < n; i++)
        {
            scal(n, sigma[i], VT + i, n);
        }
    }
    free(sigma);
    free(superb);
    return;
}

void reducedorthofactors(const double &tol, const int & m, const int &n, double * &a, double * &U, double * &VT, int &newrank)
{
    double *sigma, *superb;
    if (m <= n)
    {
        U = (double *) malloc(m * m * sizeof(double));
        VT = (double *) malloc(m * n * sizeof(double));
        sigma = (double *) malloc(m * sizeof(double));
        superb = (double *) malloc((m - 1) * sizeof(double));
        gesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, a, m, sigma, U, m, VT, m, superb);
        double s = 0;
        #pragma omp parallel for reduction (+:s)
        for (int i = 0; i < min(m,n); i++)
        {
            s += sigma[i]*sigma[i];
        }
        #pragma omp barrier
        double r = 0;
        newrank = 0;
        for (int i = m; i > 0; i--)
        {
            r += sigma[i-1]*sigma[i-1];
            if (s*tol*tol <= r * (1e0l + tol*tol))
            {
                newrank = i;
                break;
            }
        }
        for (int i = 0; i < newrank; i++)
        {
            scal(m, sigma[i], U + i * m, 1);
        }
        U = (double *) realloc(U, m * newrank * sizeof(double));
    }
    else
    {
        U = (double *) malloc(m * n * sizeof(double));
        VT = (double *) malloc(n * n * sizeof(double));
        sigma = (double *) malloc(n * sizeof(double));
        superb = (double *) malloc((n - 1) * sizeof(double));
        gesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, a, m, sigma, U, m, VT, m, superb);
        double s = 0;
        #pragma omp parallel for reduction (+:s)
        for (int i = 0; i < min(m,n); i++)
        {
            s += sigma[i]*sigma[i];
        }
        #pragma omp barrier
        double r = 0;
        newrank = 0;
        for (int i = n; i > 0; i--)
        {
            r += sigma[i-1]*sigma[i-1];
            if (s*tol*tol <= r * (1e0l + tol*tol))
            {
                newrank = i;
                break;
            }
        }
        double *tmp = (double *) malloc(newrank * n * sizeof(double));
        for (int i = 0; i < newrank; i++)
        {
            scal(n, sigma[i], VT + i, n);
            copy(n, VT + i, n, tmp, newrank);
        }
        free(VT);
        VT = tmp;
    }
    free(sigma);
    free(superb);
    return;
}

void myreducedorthofactors(const double &tol, const int & n, const int &m, double * a, double * U, double * VT, int &newrank, const int maxrank)
{
    double *sigma, *superb;

    sigma = (double *) malloc(min(m,n) * sizeof(double));
    superb = (double *) malloc((min(m,n) - 1) * sizeof(double));
    double *tmp = (double *) malloc(min(m,n) * m * sizeof(double));

    gesvd(LAPACK_COL_MAJOR, 'S', 'S', n, m, a, n, sigma, U, n, tmp, min(m,n), superb);

    long /*long*/ double s = 0;
// #pragma omp parallel for reduction (+:s)
    for (int i = min(m,n)-1; i >=0 ; i--)
    {
        s += ((long /*long*/ double)sigma[i])*((long /*long*/ double)sigma[i]);
    }
// #pragma omp barrier
    long /*long*/ double r = 0;
    newrank = 0;

    for (int i = min(n,m); i > 0; i--)
    {
        r += ((long /*long*/ double) sigma[i-1])*((long /*long*/ double) sigma[i-1]);
        if (((long /*long*/ double) s)*((long /*long*/ double) tol) * ((long /*long*/ double) tol) <= r)
        {
            newrank = i;
            break;
        }
    }
    if ((maxrank != 0) && (newrank > maxrank)) {
        newrank = maxrank;
    }

    for (int i = 0; i < newrank; i++)
    {
        scal(n, sigma[i], U + i * n, 1);
    }
    
   for (int i = 0; i < newrank; i++)
    {
        copy(m,tmp + i, min(m,n), VT + i, newrank);
    }

    free(sigma);
    free(superb);
    free(tmp);
    return;
}

int gesvd(const int &matrix_order, const char &jobu, const char &jobvt, const int &m, const int &n, double* &a, const int &lda, double* &s, double* &u, const int &ldu, double* &vt, const int &ldvt, double* &superb)
{
    return LAPACKE_dgesvd(matrix_order, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}

void copy(const int &n, const double * x, const int &incx, double *y, const int &incy)
{
    cblas_dcopy(n, x, incx, y, incy);
    return;
}

void scal(const int & n, const double & factor, double * x, const int & incx){
   cblas_dscal(n, factor, x, incx);
   return;
}

void axpy(const int &n, const double &alpha, const double *x, const int &incx, double *y, const int &incy)
{
    cblas_daxpy(n, alpha, x, incx, y, incy);
}

double nrm2(const int &n, const double *x, const int &incx)
{
    return cblas_dnrm2(n, x, incx);
}

void gemv(const CBLAS_ORDER order, const CBLAS_TRANSPOSE TransA, const int M, const int N, const double alpha, const double *A, const int lda, const double *X, const int incX, const double beta, double *Y, const int incY)
{
    cblas_dgemv(order, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
    return;
}

void gemm(const CBLAS_ORDER Order, const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, const int M, const int N, const int K, const double alpha, const double *A, const int lda, const double *B, const int ldb, const double beta, double *C, const int ldc)
{
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    return;
}

double dot_u(const int & n, const double * x, const int & incx, const double * y, const int & incy){
   return cblas_ddot(n, x, incx, y, incy);
}

double dot_c(const int & n, const double * x, const int & incx, const double * y, const int & incy){
   return cblas_ddot(n, x, incx, y, incy);
}

using namespace std;
double * pseudoinverse(double tolerance, int m, int n, double * A)
{
    if (min(m,n) == 0) return NULL;
    double * inv_A = (double *) malloc(n*m*sizeof(double));
    double * U = (double *) malloc(m*min(m, n)*sizeof(double));
    double * VT = (double *) malloc(min(m, n)*n*sizeof(double));
    double * sigma = (double *) malloc(min(m, n)*sizeof(double));
    double * superb = (double *) malloc((min(m,n)-1)*sizeof(double));
    LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, A, m, sigma, U, m, VT, min(m, n), superb);
	free(superb);
    double s = 0;
    #pragma omp parallel for reduction (+:s)
    for (int i = 0; i < min(m,n); i++)
    {
        s += sigma[i]*sigma[i];
    }
    #pragma omp barrier
    double r = 0;
    int k = 0;
    for (int i = min(m, n); i > 0; i--)
    {
        r += sigma[i-1]*sigma[i-1];
        if (s*tolerance*tolerance <= r * (1e0l + tolerance*tolerance))
        {
            k = i;
            break;
        }
    }
    #pragma omp parallel for
    for (int i = 0; i < k; i++)
    {
        cblas_dscal(n, 1e0l/sigma[i], VT + i, min(m, n));
    }
    #pragma omp barrier
    cblas_dgemm(CblasColMajor, CblasTrans, CblasTrans, n, m, k, 1e0l, VT, min(m, n), U, m, 0e0l, inv_A, n);
    free(U);
    free(VT);
    free(sigma);
    return inv_A;
}

double imag(const double & a)
{
    return 0;
}

double real(const double & a)
{
    return a;
}

int iamax(const int N, const double *X, const int incX)
{
    return cblas_idamax(N, X, incX);
}

int gesv( int matrix_order, int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb )
{
    return LAPACKE_dgesv(matrix_order, n, nrhs, a, lda, ipiv, b, ldb);
}



void ger (const CBLAS_ORDER order, const int M, const int N, const double alpha, const double *X, const int incX, const double *Y, const int incY, double *A, const int lda)
{
    cblas_dger(order, M, N, alpha, X, incX, Y, incY, A, lda);
}

void transpose(const double * mat, const int rows_number, const int column_number,  double * new_mat,  int lda_old, int lda_new)
{
    if (lda_old == 0) {
        lda_old = rows_number;        
    } else {
        if (lda_old < rows_number) {
            cout << "Warning! transpose got wrong lda! throw 1" << endl;
            throw(1);
        }            
    }
    if (lda_new == 0) {
        lda_new = column_number;        
    } else {
        if (lda_new < column_number) {
            cout << "Warning! transpose got wrong lda! throw 1" << endl;
            throw(1);
        }            
    }
    
    if (rows_number < column_number){
        double * ident = (double *) calloc (rows_number*rows_number, sizeof(double));
        double one = 1.0;
        copy(rows_number, &one, 0, ident, (rows_number + 1));
        
        gemm(CblasColMajor, CblasConjTrans, CblasNoTrans, column_number, rows_number, rows_number,
             1.0, mat,  lda_old, ident, rows_number, 0.0, new_mat, lda_new);   

        free(ident);
        // B = A^T * I
    } else {

        double * ident = (double *) calloc (column_number * column_number, sizeof(double));
        double one = 1.0;
        copy(column_number, &one, 0, ident, (column_number + 1));
        gemm(CblasColMajor, CblasNoTrans, CblasConjTrans, column_number, rows_number, column_number,
             1.0, ident, column_number, mat, lda_old, 0.0, new_mat, column_number);
        free(ident);
        // B = I * A^T
    }
}

void my_qr (int m, int n, double * a, int lda,  double * matrR, int lda_r )
{
    double* tau = (double *) malloc (min(m,n) * sizeof(double));
    LAPACKE_dgeqrf(CblasColMajor, m,  n, a, lda, tau );
    
    // copy upper triangular part  - matrix R:
    if (matrR != NULL) {
        double zero = 0.0;
        for (int i = 0; i < n; i++) {
            copy(min((i+1), m), a + i * lda, 1, matrR + i*lda_r, 1);
            if (i < m) {
              copy(m - (i + 1), &zero, 0, matrR + i * lda_r + i + 1, 1);
            }
        }
    }
    
    // genetate matrix Q:

    LAPACKE_dorgqr(CblasColMajor, m, min(m,n), min(m,n), a, lda, tau );
    
    free(tau);
    return;
}