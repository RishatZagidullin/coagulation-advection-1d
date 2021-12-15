#ifdef CUDA_FFT

#include "parallel_cross_omp.h"
#include <cassert>
using namespace convolution_helpers;

void TCross_Parallel_v1::fill_cuda_data()
{
	checkCudaErrors(cudaMalloc((void **) &d_U, this->get_rows_number() * this->rank*sizeof(double)));
	checkCudaErrors(cudaMalloc((void **) &d_V, this->get_columns_number() * this->rank*sizeof(double)));
		
	checkCudaErrors(cudaMemcpy(d_U, U, this->get_rows_number() * this->rank * sizeof(double), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_V, V, this->get_columns_number() * this->rank * sizeof(double), cudaMemcpyHostToDevice));	
}

void TCross_Parallel_v1::dealloc_cuda_data()
{
	cudaFree(d_V);
	cudaFree(d_U);    
}

void TCross_Parallel_v1::smol_conv_trapezoids(double * &x, cufftHandle & plan, cublasHandle_t &handle, double * &result)
{
//Correct the convolution with kernel A=U * transpose(V) into trapecial-weights quadrature
	int R = this->get_rank();
	int N = this->rows_number;

	dim3 block(128);
	dim3 grid((N+block.x-1)/block.x);

	//double *d_res;
	
	this->smol_conv(x, plan, handle, result);

	cudaMemset(result, 0.0, sizeof(double));
	smol_conv_help_function<<<grid, block>>>(x, d_U, d_V, result, N, R);
	//return d_res;
}

void TCross_Parallel_v1::smol_conv(double * &x, cufftHandle & plan, cublasHandle_t &handle, double * &result)
{

// allocate memory
	int R = this->get_rank();
	int N = this->rows_number;
	int M = this->columns_number;

	dim3 block(128);
	dim3 grid((N+block.x-1)/block.x);

	if (R <= 0) {cout << "unacceptable RANK" <<endl; return;}

	double ones[R];
	for (int j = 0; j < R; j++)
	{
		ones[j] = 1.0l;
	}
	double *d_ones;
	cudaMalloc((void **) &d_ones, R * sizeof(double));
	cudaMemcpy(d_ones, ones, R*sizeof(double), cudaMemcpyHostToDevice);

       	double *temp_res;
	cudaMalloc((void **) &temp_res, R*N*sizeof(double));
       	//double *result;
	//cudaMalloc((void **) &result, N*sizeof(double));

	int mem_size = sizeof(cufftDoubleComplex) * N * R;

	cufftDoubleComplex *d_ub;
	cudaMalloc((void **)&d_ub, mem_size);

	cufftDoubleComplex *d_vb;
	cudaMalloc((void **)&d_vb, mem_size);

	for (int r = 0; r < R; r++)
      	{
		MulAndComplexCast<<<grid, block >>>(&d_ub[r*M], &d_U[r*M], x, N);
		MulAndComplexCast<<<grid, block >>>(&d_vb[r*M], &d_V[r*N], x, M);
	}
		
	cufftExecZ2Z(plan, (d_vb), (d_vb), CUFFT_FORWARD);	
	cufftExecZ2Z(plan, (d_ub), (d_ub), CUFFT_FORWARD);
		
	ComplexPointwiseMulAndScale<<<grid, block >>>(d_ub, d_vb, N * R, 1.0l / N);

	cufftExecZ2Z(plan, (d_ub), (d_ub), CUFFT_INVERSE);
		
	double alpha = 1.0;
	double beta = 0.0;
	for (int i = 0; i < R; i++)
	{
		DoubleCastFromComplex<<<grid, block>>>(&d_ub[i*N], &temp_res[i*N], N);
	}
	cublasDgemv(handle, CUBLAS_OP_N, N, R, &alpha, temp_res, N, d_ones, 1, &beta, result, 1);

	cudaFree(d_vb);
	cudaFree(d_ub);
	cudaFree(d_ones);
	cudaFree(temp_res);
}

void TCross_Parallel_v1::matvec(double* &x, cublasHandle_t &handle, double* &result, const char &option)
{
	int R = this->get_rank();
	int M = this->rows_number;
	int N = this->columns_number;

	if (option == 'f')
	{
		//full matvec
		double alpha = 1.0;
		double beta = 0.0;

		double *d_buff;
		cudaMalloc((void **) &d_buff, R*sizeof(double));
		// compute transpose(V)*x
                cublasDgemv(handle, CUBLAS_OP_T, N, R, &alpha, d_V, N, x, 1, &beta, d_buff, 1);
		// compute U * (V*x)

		cublasDgemv(handle, CUBLAS_OP_N, M, R, &alpha, d_U, M, d_buff, 1, &beta, result, 1);

		cudaFree(d_buff);
	}
	else
	{
		printf("incorrect matvec option\nreturn null\n");
	}
}
#endif
