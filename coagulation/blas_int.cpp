#include "blas.h"
void copy(const int & n, const int * x, const int & incx, int * y, const int & incy){
    if (n == 0) return;
    #pragma omp parallel for
    for (int i = 0; i < n; i++){
        int ix = i * incx;
        int iy = i * incy;
        y[iy] = x[ix];
    }
    #pragma omp barrier
    return;
}
