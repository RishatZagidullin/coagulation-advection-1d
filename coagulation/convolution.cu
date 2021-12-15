#ifdef CUDA_FFT
#include "convolution.cuh"
namespace convolution_helpers
{
	__device__ __host__ cufftDoubleComplex ComplexScale(cufftDoubleComplex a, double s)
	{
		cufftDoubleComplex c;
		c.x = s * a.x;
		c.y = s * a.y;
		return c;
	}

	__device__ __host__ cufftDoubleComplex ComplexMul(cufftDoubleComplex a, cufftDoubleComplex b)
	{
		cufftDoubleComplex c;
		c.x = a.x * b.x - a.y * b.y;
		c.y = a.x * b.y + a.y * b.x;
		return c;
	}

	__global__ void ComplexPointwiseMulAndScale(cufftDoubleComplex *a, const cufftDoubleComplex *b, int size, double scale)
	{
		const int numThreads = blockDim.x * gridDim.x;
		const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
		for (int i = threadID; i < size; i += numThreads)
		{
			a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);
		}
	}

	__global__ void MulAndComplexCast(cufftDoubleComplex *result, const double *x, const double *coef, int size)
	{
		const int numThreads = blockDim.x * gridDim.x;
		const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
		for (int i = threadID; i < size; i += numThreads)
		{
			result[i].x = coef[i] * x[i];
			result[i].y = 0.0;
		}
	}

	__global__ void DoubleCastFromComplex(const cufftDoubleComplex *x, double *result, int size)
	{
		const int numThreads = blockDim.x * gridDim.x;
		const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
		for (int i = threadID; i < size; i += numThreads)
		{
			result[i] = x[i].x;
		}
	}

	__global__ void smol_conv_help_function(const double *x, const double *U, const double *V, double *result, const int N, const int R)
	{
		const int numThreads = blockDim.x * gridDim.x;
		const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
		for (int i = threadID; i < N; i += numThreads)
		{
			if(i!=0) for(int j = 0; j < R; j++) result[i] -= 0.5 * (U[N * j + i] * V[N * j] * x[0] * x[i] + U[N * j] * V[N * j + i] * x[i] * x[0]);
		}
	}
}
#endif
