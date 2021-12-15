#pragma once

#ifdef CUDA_FFT
#include <cuda_runtime.h>
#include <cufft.h>
#include <cublas_v2.h>
#include "../cuda_utils/helper_cuda.h"

namespace convolution_helpers
{
	__device__ __host__ cufftDoubleComplex ComplexScale(cufftDoubleComplex, double);
	__device__ __host__ cufftDoubleComplex ComplexMul(cufftDoubleComplex, cufftDoubleComplex);
	__global__ void ComplexPointwiseMulAndScale(cufftDoubleComplex *, const cufftDoubleComplex *, int, double);
	__global__ void MulAndComplexCast(cufftDoubleComplex *result, const double *x, const double *coef, int size);
	__global__ void DoubleCastFromComplex(const cufftDoubleComplex *x, double *result, int size);
	__global__ void smol_conv_help_function(const double *x, const double *U, const double *V, double *result, const int N, const int R);
}
#endif
