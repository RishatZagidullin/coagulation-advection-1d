#ifndef PARALLEL_CROSS_H
#define PARALLEL_CROSS_H
#include "cross.h"
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "blas.h"
#include <omp.h>
#include "maxvol.h"
#include "stdio.h"

//**************************RISHAT*ADDED*BEGIN*************************
#ifdef FFTW
#include "convolution.h"
#endif

#ifdef CUDA_FFT
#include "convolution.cuh"
#endif

#ifdef MKL_FFT
#include "mkl_dfti.h"
#endif
//**************************RISAT*ADDED*END***************************

struct TCross_Parallel_v1_Parameters
{
    double tolerance;
    int maximal_iterations_number, number_of_checked_elements, max_rank, rank_increase, stop_rank;
    bool memory_strategy, start_from_column;
    TCross_Parallel_v1_Parameters();
};

struct TCross_Parallel_v1_Work_Data{
    bool *J, *K;
    double *current_row, *current_column;
    int omp_column_threads_num, omp_row_threads_num, *omp_column_start, *omp_column_num, *omp_row_start, *omp_row_num;
    TVolume max_volume;
    double global_max;
    TMatrix * matrix;
    TDifferenceMatrix work_matrix;
    TCross_Parallel_v1_Parameters parameters;
    TCross_Parallel_v1_Work_Data(TMatrix *, TMatrix *);
    ~TCross_Parallel_v1_Work_Data();
};

class TCross_Parallel_v1: public TCross_Base<TCross_Parallel_v1_Work_Data, TCross_Parallel_v1_Parameters>
{
    private:
        double *U, *V, *C, *hat_A_inv, *RT, tolerance, norm, *AR, *CAT;
	//**************************RISHAT*ADDED*BEGIN*************************
#ifdef CUDA_FFT
	double *d_U, *d_V;
	void fill_cuda_data();
	void dealloc_cuda_data();
#endif
	//**************************RISHAT*ADDED*END***************************
        int * rows_numbers, * columns_numbers;
        void Prepare_Data(TCross_Parallel_v1_Work_Data &, const TCross_Parallel_v1_Parameters &);
        void Search_Max_Volume(TCross_Parallel_v1_Work_Data &);
        bool Stopping_Criteria(TCross_Parallel_v1_Work_Data &);
        void Update_Cross(TCross_Parallel_v1_Work_Data &);
        void get_diff_column(const int &, const TCross_Parallel_v1_Work_Data &, double *&);
        void get_diff_row(const int &, const TCross_Parallel_v1_Work_Data &, double *&);
	//**************************RISHAT*ADDED*BEGIN*************************
#ifdef CUDA_FFT
	void smol_conv(double* &, cufftHandle &, cublasHandle_t &, double * &);
#endif
#ifdef FFTW
	double *smol_conv(double* &, fftw_complex * & ub, fftw_complex * & uv, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse);
#endif
#ifdef MKL_FFT
	double *smol_conv(double* &, VSLConvTaskPtr *&);
#endif
	//**************************RISHAT*ADDED*END***************************
    public:
        void get_approximations(double *, double *);

	int get_columns_amount();
	//**************************RISHAT*ADDED*BEGIN*************************
#ifdef CUDA_FFT
        void matvec(double *&x, cublasHandle_t &, double * &, const char &option='f');
        void smol_conv_trapezoids(double *&, cufftHandle &, cublasHandle_t &, double *&);
#endif
#ifndef CUDA_FFT
	double *matvec(double *&x, const char &option='f');
#endif
#ifdef FFTW
        double *smol_conv_trapezoids(double *&, fftw_complex * & ub, fftw_complex * & uv, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse);
#endif
#ifdef MKL_FFT
	double *smol_conv_trapezoids(double *&, VSLConvTaskPtr *&);
#endif
	//**************************RISHAT*ADDED*END*************************
        double *smol_conv_discrete(double *&, fftw_complex * & ub, fftw_complex * & uv, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse);
        double value(const int &, const int &);
        int get_row_number(const int &) const;
        int get_column_number(const int &) const;
        TCross_Parallel_v1();
        ~TCross_Parallel_v1();
        const double *export_C();
        const double *export_hat_A_inv();
        const double *export_RT();
        const double *export_CAT();
        const double *export_AR();
};

#endif
