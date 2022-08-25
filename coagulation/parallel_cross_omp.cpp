#include "parallel_cross_omp.h"
#include <cassert>
using namespace convolution_helpers;
bool TCross_Parallel_v1::Stopping_Criteria(TCross_Parallel_v1_Work_Data &work_data){
    if (abs(work_data.max_volume.get_volume()) > abs(work_data.global_max))
    {
        work_data.global_max = abs(work_data.max_volume.get_volume());
    }
    if (((work_data.parameters.stop_rank > 0) && (this->rank >= work_data.parameters.stop_rank)) ||(sqrt(real(norm))*abs(tolerance) >= abs(work_data.max_volume.get_volume())*sqrt(this->get_columns_number()-this->rank)*sqrt(this->get_rows_number()-this->rank)))
    {
        #pragma omp parallel sections
        {
            #pragma omp section
                U = (double *) realloc(U, this->rank * this->get_rows_number() * sizeof(double));
            #pragma omp section
                V = (double *) realloc(V, this->rank * this->get_columns_number() * sizeof(double));
            #pragma omp section
                C = (double *) malloc(this->rank * this->get_rows_number() * sizeof(double));
            #pragma omp section
                RT = (double *) malloc(this->rank * this->get_columns_number() * sizeof(double));
            #pragma omp section
                rows_numbers = (int *) realloc(rows_numbers ,this->rank*sizeof(int));
            #pragma omp section
                columns_numbers = (int *) realloc(columns_numbers,this->rank*sizeof(int));
            #pragma omp section
                CAT = (double *) malloc(this->rank * this->get_rows_number() * sizeof(double));
            #pragma omp section
                AR = (double *) malloc(this->rank * this->get_columns_number() * sizeof(double));
        }
        #pragma omp barrier
        maxvol(rows_number, rank, U, rows_numbers, tolerance, CAT);
        for (int k = 0; k < rank; k++)
        {
            gemv(CblasColMajor, CblasNoTrans, get_columns_number(), rank, 1.0, V, get_columns_number(), U + rows_numbers[k], get_rows_number(), 0.0, RT + k * get_columns_number(), 1);
        }
        maxvol(columns_number, rank, V, columns_numbers, tolerance, AR);
        for (int k = 0; k < rank; k++)
        {
            gemv(CblasColMajor, CblasNoTrans, get_rows_number(), rank, 1.0, U, get_rows_number(), V + columns_numbers[k], get_columns_number(), 0.0, C + k * get_rows_number(), 1);
        }

	//*******************************RISHAT*ADDED*BEGIN****************************************
#ifdef CUDA_FFT
	fill_cuda_data();
#endif
	//******************************RISHAT*ADDED*END*******************************************
        int *tmp_ipiv;
        double *tmp_matrix;
        tmp_matrix = (double *)malloc(rank * rank * sizeof(double));
        tmp_ipiv = (int *)malloc(rank * sizeof(int));
        for (int i = 0; i < rank; i++)
        {
                copy(rows_number, C + i*rows_number, 1, CAT + i, rank);
        }
        for (int i = 0; i < rank ; i++)
        {
                double *pointer1 = C + rows_numbers[i];
                double *pointer2 = tmp_matrix+i*rank;
                copy(rank, pointer1, rows_number, pointer2, 1);
        }
        gesv(LAPACK_COL_MAJOR, rank, rows_number, tmp_matrix, rank, tmp_ipiv, CAT, rank);
        for (int i = 0; i < rank; i++)
        {
                copy(columns_number, RT + i*columns_number, 1, AR + i, rank);
        }
        for (int i = 0; i < rank ; i++)
        {
                double *pointer1 = RT + columns_numbers[i];
                double *pointer2 = tmp_matrix + i*rank;
                copy(rank, pointer1, columns_number, pointer2, 1);
        }
        gesv(LAPACK_COL_MAJOR, rank, columns_number, tmp_matrix, rank, tmp_ipiv, AR, rank);
        free(tmp_ipiv);
        free(tmp_matrix);
        double *hat_A = (double *) malloc(this->rank * this->rank * sizeof(double));
        #pragma omp parallel for
        for (int k = ((int) 0); k < this->rank; k++)
        {
            copy(this->rank, RT + columns_numbers[k], this->get_columns_number(), hat_A + k*this->rank, 1);
        }
        hat_A_inv = pseudoinverse(1e-8, this->rank, this->rank, hat_A);
        free(hat_A);
        work_data.current_row = NULL;
        work_data.current_column = NULL;
        return true;
    }
    else
    {
        return false;
    }
}


TCross_Parallel_v1_Parameters::TCross_Parallel_v1_Parameters():tolerance(0.0), maximal_iterations_number(1), number_of_checked_elements(1), max_rank(1), rank_increase(1), memory_strategy(true), stop_rank(0), start_from_column(true)
{
}

const double * TCross_Parallel_v1::export_C()
{
    return C;
}

const double * TCross_Parallel_v1::export_CAT()
{
    return CAT;
}

const double * TCross_Parallel_v1::export_AR()
{
    return AR;
}

const double * TCross_Parallel_v1::export_hat_A_inv()
{
    return hat_A_inv;
}

const double * TCross_Parallel_v1::export_RT()
{
    return RT;
}

using namespace std;

void TCross_Parallel_v1::get_diff_column(const int &j, const TCross_Parallel_v1_Work_Data & work_data, double *&result)
{
    int i;
    #pragma omp parallel for
    for (i = ((int) 0); i < this->rows_number; i++)
    {
        if (!work_data.K[i])
        {
            result[i] = work_data.matrix->value(i, j);
        }
    }
    gemv(CblasColMajor, CblasNoTrans, this->get_rows_number(), this->rank, - 1.0, U, this->get_rows_number(), V + j, this->get_columns_number(), 1.0, result, 1);
    return;
}

void TCross_Parallel_v1::get_diff_row(const int &i, const TCross_Parallel_v1_Work_Data & work_data, double *&result)
{
    int j;
    #pragma omp parallel for
    for (j = ((int) 0); j < this->columns_number; j++)
    {
        if (!work_data.J[j])
        {
            result[j] = work_data.matrix->value(i, j);
        }
    }
    gemv(CblasColMajor, CblasNoTrans, this->get_columns_number(), this->rank, - 1.0, V, this->get_columns_number(), U + i, this->get_rows_number(), 1.0, result, 1);
    return;
}

TCross_Parallel_v1_Work_Data::TCross_Parallel_v1_Work_Data(TMatrix *original_matrix, TMatrix *approximated_matrix): max_volume(1), work_matrix(original_matrix, approximated_matrix)
{
    K = (bool *) malloc(original_matrix->get_rows_number()*sizeof(bool));
    #pragma omp parallel for
    for (int i = ((int) 0); i < original_matrix->get_rows_number(); i++) K[i] = false;
    J = (bool *) malloc(original_matrix->get_columns_number()*sizeof(bool));
    #pragma omp parallel for
    for (int j = ((int) 0); j < original_matrix->get_columns_number(); j++) J[j] = false;
    global_max = 0.0;
    matrix = original_matrix;
    omp_column_threads_num = omp_get_max_threads();
    omp_row_threads_num = omp_column_threads_num;
    int avg_column = (original_matrix->get_rows_number() - 1) / omp_column_threads_num + 1;
    int avg_row = (original_matrix->get_columns_number() - 1) / omp_row_threads_num + 1;
    omp_column_threads_num = (original_matrix->get_rows_number() - 1) / avg_column + 1;
    omp_row_threads_num = (original_matrix->get_columns_number() - 1) / avg_row + 1;
    omp_column_start = (int *) malloc(omp_column_threads_num * sizeof(int));
    omp_column_num = (int *) malloc(omp_column_threads_num * sizeof(int));
    omp_row_start = (int *) malloc(omp_row_threads_num * sizeof(int));
    omp_row_num = (int *) malloc(omp_row_threads_num * sizeof(int));
    omp_column_start[((int) 0)] = ((int) 0);
    omp_row_start[((int) 0)] = ((int) 0);
    for (int i = 1; i < omp_column_threads_num; i++)
    {
        omp_column_start[i] = omp_column_start[i - 1] + avg_column;
        omp_column_num[i - 1] = avg_column;
    }
    omp_column_num[omp_column_threads_num - 1] = original_matrix->get_rows_number() - omp_column_start[omp_column_threads_num - 1];
    for (int i = 1; i < omp_row_threads_num; i++)
    {
        omp_row_start[i] = omp_row_start[i - 1] + avg_row;
        omp_row_num[i - 1] = avg_row;
    }
    omp_row_num[omp_row_threads_num - 1] = original_matrix->get_columns_number() - omp_row_start[omp_row_threads_num - 1];
}

TCross_Parallel_v1_Work_Data::~TCross_Parallel_v1_Work_Data()
{
    free(K);
    free(J);
    free(current_row);
    free(current_column);
    free(omp_column_start);
    free(omp_column_num);
    free(omp_row_start);
    free(omp_row_num);
}

void TCross_Parallel_v1::Prepare_Data(TCross_Parallel_v1_Work_Data &work_data,  const TCross_Parallel_v1_Parameters & parameters)
{
    norm = 0.0;
    free(U);
    free(V);
    free(C);
    free(hat_A_inv);
    free(RT);
    free(CAT);
    free(AR);
    U = NULL;
    V = NULL;
    C = NULL;
    CAT = NULL;
    AR = NULL;
    hat_A_inv = NULL;
    RT = NULL;
    free(rows_numbers);
    free(columns_numbers);
    rows_numbers = NULL;
    columns_numbers = NULL;
    tolerance = parameters.tolerance;
    work_data.parameters = parameters;
    if (work_data.parameters.max_rank <= ((int) 0))
    {
        work_data.parameters.max_rank = 1;
    }
    if ((work_data.parameters.rank_increase) <= ((int) 0) && (work_data.parameters.memory_strategy))
    {
        work_data.parameters.rank_increase = 1;
    }
    else if ((work_data.parameters.rank_increase) <= 1 && (!work_data.parameters.memory_strategy))
    {
        work_data.parameters.rank_increase = 1 + 1;
    }
    #pragma omp parallel sections
    {
        #pragma omp section
            U = (double *) malloc(work_data.parameters.max_rank * this->get_rows_number() * sizeof(double));
        #pragma omp section
            V = (double *) malloc(work_data.parameters.max_rank * this->get_columns_number() * sizeof(double));
        #pragma omp section
            rows_numbers = (int *) realloc(rows_numbers, work_data.parameters.max_rank*sizeof(int));
        #pragma omp section
            columns_numbers = (int *) realloc(columns_numbers, work_data.parameters.max_rank*sizeof(int));
    }
    this->rank = ((int) 0);
}

TCross_Parallel_v1::TCross_Parallel_v1(): TCross_Base<TCross_Parallel_v1_Work_Data, TCross_Parallel_v1_Parameters>()
{
    U = NULL;
    V = NULL;
    C = NULL;
    hat_A_inv = NULL;
    RT = NULL;
    CAT = NULL;
    AR = NULL;
    rows_numbers = NULL;
    columns_numbers = NULL;
    this->rank = ((int) 0);
};

double TCross_Parallel_v1::value(const int &i, const int &j)
{
    return dot_u(this->rank, U + i, this->get_rows_number(), V + j, this->get_columns_number());
}

void TCross_Parallel_v1::Search_Max_Volume(TCross_Parallel_v1_Work_Data &work_data)
{
    if (this->rank == min(this->get_rows_number(), this->get_columns_number())) return;
    if (this->rank >= work_data.parameters.max_rank)
    {
        if (work_data.parameters.memory_strategy)
        {
            work_data.parameters.max_rank += work_data.parameters.rank_increase;
        }
        else
        {
            work_data.parameters.max_rank *= work_data.parameters.rank_increase;
        }
        #pragma omp parallel sections
        {
            #pragma omp section
                U = (double *) realloc(U, work_data.parameters.max_rank * this->get_rows_number() * sizeof(double));
            #pragma omp section
                V = (double *) realloc(V, work_data.parameters.max_rank * this->get_columns_number() * sizeof(double));
            #pragma omp section
                rows_numbers = (int *) realloc(rows_numbers, work_data.parameters.max_rank*sizeof(int));
            #pragma omp section
                columns_numbers = (int *) realloc(columns_numbers, work_data.parameters.max_rank*sizeof(int));
        }
        #pragma omp barrier
    }
    work_data.current_column = U + this->rank * this->get_rows_number();
    work_data.current_row = V + this->rank * this->get_columns_number();
    int ik = rand() % (this->get_rows_number() - this->rank), jk = rand() % (this->get_columns_number() - this->rank);
    int i(0), j(0);
    int ia(0), ja(0);
    #pragma omp parallel sections
    {
        #pragma omp section
            while (ia <= ik)
            {
                if (!work_data.K[i]) ia++;
                i++;
            }
        #pragma omp section
            while (ja <= jk)
            {
                if (!work_data.J[j]) ja++;
                j++;
            }
    }
    i--;
    j--;
    #pragma omp barrier
    work_data.max_volume.set(work_data.work_matrix.value(i, j), &i, &j);
    int iters = ((int) 0), prev_column(j), prev_row(i);
    bool b;
    if (work_data.parameters.start_from_column)
    {
        do
        {
            b = (work_data.max_volume.get_column_position() == prev_column);
            if (iters == 0)
            {
                b = false;
            }
            if (!b)
            {
                j = work_data.max_volume.get_column_position();
                get_diff_column(j, work_data, work_data.current_column);
                i = work_data.max_volume.get_row_position();
                int ind[work_data.omp_column_threads_num];
                #pragma omp parallel for
                for (int o = ((int) 0); o < work_data.omp_column_threads_num; o++)
                {
                    ind[o] = i;
                    for (int k = work_data.omp_column_start[o]; k < work_data.omp_column_start[o] + work_data.omp_column_num[o]; k++)
                    {
                        if (!work_data.K[k])
                        {
                            if (abs(work_data.current_column[ind[o]]) < abs(work_data.current_column[k]))
                            {
                                ind[o] = k;
                            }
                        }
                    }
                }
                #pragma omp barrier
                i = ind[((int) 0)];
                for (int o = 1; o < work_data.omp_column_threads_num; o++)
                {
                    if (abs(work_data.current_column[i]) < abs(work_data.current_column[ind[o]]))
                    {
                        i = ind[o];
                    }
                }
                work_data.max_volume.set(work_data.current_column[i], &i, &j);
                prev_column = work_data.max_volume.get_column_position();
            }
            b = (work_data.max_volume.get_row_position() == prev_row);
            if (iters == 0)
            {
                b = false;
            }
            if (!b)
            {
                i = work_data.max_volume.get_row_position();
                get_diff_row(i, work_data, work_data.current_row);
                j = work_data.max_volume.get_column_position();
                int ind[work_data.omp_row_threads_num];
                #pragma omp parallel for
                for (int o = ((int) 0); o < work_data.omp_row_threads_num; o++)
                {
                    ind[o] = j;
                    for (int k = work_data.omp_row_start[o]; k < work_data.omp_row_start[o] + work_data.omp_row_num[o]; k++)
                    {
                        if (!work_data.J[k])
                        {
                            if (abs(work_data.current_row[ind[o]]) < abs(work_data.current_row[k]))
                            {
                                ind[o] = k;
                            }
                        }
                    }
                }
                #pragma omp barrier
                j = ind[((int) 0)];
                for (int o = 1; o < work_data.omp_row_threads_num; o++)
                {
                    if (abs(work_data.current_row[j]) < abs(work_data.current_row[ind[o]]))
                    {
                        j = ind[o];
                    }
                }
                work_data.max_volume.set(work_data.current_row[j], &i, &j);
                prev_row = work_data.max_volume.get_row_position();
            }
            iters++;
        } while ((!b) || ((work_data.parameters.maximal_iterations_number > 0) && (iters < work_data.parameters.maximal_iterations_number)));
    }
    else
    {
        do
        {
            b = (work_data.max_volume.get_row_position() == prev_row);
            if (iters == 0)
            {
                b = false;
            }
            if (!b)
            {
                i = work_data.max_volume.get_row_position();
                get_diff_row(i, work_data, work_data.current_row);
                j = work_data.max_volume.get_column_position();
                int ind[work_data.omp_row_threads_num];
                #pragma omp parallel for
                for (int o = ((int) 0); o < work_data.omp_row_threads_num; o++)
                {
                    ind[o] = j;
                    for (int k = work_data.omp_row_start[o]; k < work_data.omp_row_start[o] + work_data.omp_row_num[o]; k++)
                    {
                        if (!work_data.J[k])
                        {
                            if (abs(work_data.current_row[ind[o]]) < abs(work_data.current_row[k]))
                            {
                                ind[o] = k;
                            }
                        }
                    }
                }
                #pragma omp barrier
                j = ind[((int) 0)];
                for (int o = 1; o < work_data.omp_row_threads_num; o++)
                {
                    if (abs(work_data.current_row[j]) < abs(work_data.current_row[ind[o]]))
                    {
                        j = ind[o];
                    }
                }
                work_data.max_volume.set(work_data.current_row[j], &i, &j);
                prev_row = work_data.max_volume.get_row_position();
            }
            b = (work_data.max_volume.get_column_position() == prev_column);
            if (iters == 0)
            {
                b = false;
            }
            if (!b)
            {
                j = work_data.max_volume.get_column_position();
                get_diff_column(j, work_data, work_data.current_column);
                i = work_data.max_volume.get_row_position();
                int ind[work_data.omp_column_threads_num];
                #pragma omp parallel for
                for (int o = ((int) 0); o < work_data.omp_column_threads_num; o++)
                {
                    ind[o] = i;
                    for (int k = work_data.omp_column_start[o]; k < work_data.omp_column_start[o] + work_data.omp_column_num[o]; k++)
                    {
                        if (!work_data.K[k])
                        {
                            if (abs(work_data.current_column[ind[o]]) < abs(work_data.current_column[k]))
                            {
                                ind[o] = k;
                            }
                        }
                    }
                }
                #pragma omp barrier
                i = ind[((int) 0)];
                for (int o = 1; o < work_data.omp_column_threads_num; o++)
                {
                    if (abs(work_data.current_column[i]) < abs(work_data.current_column[ind[o]]))
                    {
                        i = ind[o];
                    }
                }
                work_data.max_volume.set(work_data.current_column[i], &i, &j);
                prev_column = work_data.max_volume.get_column_position();
            }
            iters++;
        } while ((!b) || ((work_data.parameters.maximal_iterations_number > 0) && (iters < work_data.parameters.maximal_iterations_number)));
    }
    if (prev_column != work_data.max_volume.get_column_position())
    {
        get_diff_column(work_data.max_volume.get_column_position(), work_data, work_data.current_column);
    }
    if (prev_row != work_data.max_volume.get_row_position())
    {
        get_diff_row(work_data.max_volume.get_row_position(), work_data, work_data.current_row);
    }
}

void TCross_Parallel_v1::Update_Cross(TCross_Parallel_v1_Work_Data &work_data)
{
    double column_factor(1.0/sqrt(abs(work_data.max_volume.get_volume()))), row_factor = 1.0 / (work_data.max_volume.get_volume() * column_factor);
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            #pragma omp parallel for
            for (int i = ((int) 0); i < this->get_rows_number(); i++)
            {
                if (!work_data.K[i])
                {
                    work_data.current_column[i] *= column_factor;
                }
                else
                {
                    work_data.current_column[i] = 0.0;
                }
            }
        }
        #pragma omp section
        {
            #pragma omp parallel for
            for (int j = ((int) 0); j < this->get_columns_number(); j++)
            {
                if (!work_data.J[j])
                {
                    work_data.current_row[j] *= row_factor;
                }
                else
                {
                    work_data.current_row[j] = 0.0;
                }
            }
        }
    }
    work_data.K[work_data.max_volume.get_row_position()] = true;
    work_data.J[work_data.max_volume.get_column_position()] = true;
    rows_numbers[this->rank] = work_data.max_volume.get_row_position();
    columns_numbers[this->rank] = work_data.max_volume.get_column_position();
    double *b = (double *) malloc ((this->rank + 1) * sizeof(double)), *x = (double *) malloc ((this->rank + 1) * sizeof(double));
    gemv(CblasColMajor, CblasConjTrans, this->get_rows_number(), this->rank + 1, 1.0, U, this->get_rows_number(), work_data.current_column, 1, 0.0, b, 1);
    gemv(CblasColMajor, CblasConjTrans, this->get_columns_number(), this->rank + 1, 1.0, V, this->get_columns_number(), work_data.current_row, 1, 0.0, x, 1);
    norm += 2.0*real(dot_c(this->rank, b, 1, x, 1)) + b[this->rank]*x[this->rank];
    free(x);
    free(b);
    this->rank++;
    work_data.current_column = NULL;
    work_data.current_row = NULL;
}

int TCross_Parallel_v1::get_row_number(const int & k) const
{
    return rows_numbers[k];
}

int TCross_Parallel_v1::get_column_number(const int & k) const
{
    return columns_numbers[k];
}

//*******************************RISHAT*ADDED*BEGIN****************************************
#ifdef FFTW
double *TCross_Parallel_v1::smol_conv(double * &x, fftw_complex * & ub, fftw_complex * & vb, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse)
{
	int R = this->get_rank();
	int N = this->rows_number;
	int M = this->columns_number;

	if (R <= 0) return NULL;

	double ones[R];
	for (int j = 0; j < R; j++)
	{
		ones[j] = 1.0l;
	}
       	double *tmp_res = (double *) malloc(R * N * sizeof(double));
       	double *result = (double *) malloc(N * sizeof(double));
	#pragma omp parallel for collapse(2)
      	for (int r = 0; r < R; r++)
      	{
		for (int i = 0; i < N; i++)
		{
			vb[i + r * N][0] = V[i + r * N] * x[i];
			vb[i + r * N][1] = 0.0;
		}
	}
	#pragma omp parallel for collapse(2)
	for (int r = 0; r < R; r++)
      	{
		for (int i = 0; i < M; i++)
		{
			ub[i + r * M][0] = U[i + r * M] * x[i];
			ub[i + r * M][1] = 0.0;
		}
	}
	for (int i = 0; i < R; i++)
	{
		fftw_execute(plan_v[i]);
		fftw_execute(plan_u[i]);
		
	}
	ComplexPointwiseMulAndScale_fftw(ub, vb, N * R, 1.0l / N);
	for (int i = 0; i < R; i++)
	{
		fftw_execute(plan_inverse[i]);
		
	}	
	#pragma omp parallel for collapse(2)
	for (int i = 0; i < R; i++)
	{
		for (int j = 0; j < N; j++)
		{
			tmp_res[i * N + j] = ub[i * N + j][0];
		}
	}
	gemv(CblasColMajor, CblasNoTrans, N, R, 1.0l, tmp_res, N, ones, 1, 0.0l, result, 1);
	free(tmp_res);
	return result;
}

double *TCross_Parallel_v1::smol_conv_trapezoids(double * &x, fftw_complex * & ub, fftw_complex * & uv, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse)
{
	int R = this->get_rank();
	int N = this->rows_number;
	double *result = (double *)malloc (N *sizeof(double));
	free(result);
	result = this->smol_conv(x, ub, uv, plan_v, plan_u, plan_inverse);
	result[0] = 0.0;
	#pragma omp parallel for
	for (int i = int (1); i < N; i++)
		for(int j = 0; j < R; j++)
			result[i] -= 0.5 * (U[N * j + i] * V[N * j] * x[0] * x[i] + U[N * j] * V[N * j + i] * x[i] * x[0]);
	result[0] = 0.0;
	return result;
}
#endif
#ifdef MKL_FFT
double *TCross_Parallel_v1::smol_conv(double * &x, VSLConvTaskPtr * &tasks)
{
	int R = this->get_rank();
	int N = this->rows_number;
	int M = this->columns_number;
        if (R <= 0) return NULL;

        double ones[R];
	for (int j = 0; j < R; j++)
	{
		ones[j] = 1.0l;
        }
        assert(M==N);
        int status;
        double t = omp_get_wtime();
        double *ub = (double *) malloc(R * M * sizeof(double));
        double *vb = (double *) malloc(R * N * sizeof(double));\
        double *tmp_res = (double *) malloc(R * N * sizeof(double));
        double *result = (double *) malloc(N * sizeof(double));
        #pragma omp parallel for collapse(2)
        for (int r = 0; r < R; r++)
        {
		for (int i = 0; i < N; i++)
		{
			ub[i + r * M] = U[i + r * M] * x[i];
			vb[i + r * N] = V[i + r * N] * x[i];
		}
	}
        #pragma omp parallel for
        for (int r = 0; r < R; r++)
        {
		status = vsldConvExec1D(tasks[r], ub + r * M, 1, vb + r * N, 1, tmp_res + r * N, 1);
        }
        gemv(CblasColMajor, CblasNoTrans, N, R, 1.0l, tmp_res, N, ones, 1, 0.0l, result, 1);
        t = omp_get_wtime() - t;
        free(ub);
        free(vb);
        free(tmp_res);
        return result;
}

double *TCross_Parallel_v1::smol_conv_trapezoids(double * &x, VSLConvTaskPtr * &tasks)
{
	int R = this->get_rank();
	int N = this->rows_number;
	double *result = (double *)malloc (N *sizeof(double));
	free(result);
	result = this->smol_conv(x, tasks);
	result[0] = 0.0;
	#pragma omp parallel for
	for (int i = int (1); i < N; i++)
		for(int j = 0; j < R; j++)
			result[i] -= 0.5 * (U[N * j + i] * V[N * j] * x[0] * x[i] + U[N * j] * V[N * j + i] * x[i] * x[0]);
	result[0] = 0.0;
	return result;
}
#endif
//*******************************RISHAT*ADDED*END******************************************

double *TCross_Parallel_v1::smol_conv_discrete(double *&x, fftw_complex * & ub, fftw_complex * & uv, fftw_plan * & plan_v, fftw_plan *& plan_u, fftw_plan *& plan_inverse)
{
	int N = this->rows_number;
	double *result = (double *)malloc (N * sizeof(double));
	free(result);
	result = this ->smol_conv(x, ub, uv, plan_v, plan_u, plan_inverse);
	double bubble1, bubble2;
	bubble1 = result[0];
	for (int i = 0; i < (N - 1); i++)
	{
		bubble2 = result[i + 1];
		result[i + 1] = bubble1;
		bubble1 = bubble2;
	}
	result[0] = 0;
	return result;
}
#ifndef CUDA_FFT
double* TCross_Parallel_v1::matvec(double* &x, const char &option)
{
	int R = this->get_rank();
	int M = this->rows_number;
	int N = this->columns_number;
	double *result = (double*) malloc( M * sizeof(double) );
	double w;

	if (option == 'f')
	{
		double *buff  = (double*) malloc( R * sizeof(double) );
                gemv(CblasColMajor, CblasTrans, N, R, 1.0l, V, N, x, 1, 0.0l, buff, 1);

                gemv(CblasColMajor, CblasNoTrans, M, R, 1.0l, U, M, buff, 1, 0.0l, result, 1);
		free(buff);
		buff = NULL;

	}
	else if (option == 'p')
	{
                #pragma omp parallel for
		for (int i = 0; i < M; i++)
			result[i] = 0.0;
		for (int i = 0; i < R; i++)
		{
			w = 0.0;
			for (int j = 0; j < M; j++)
			{
				if (j < N)
					w += V[N * i + (N - 1 - j)] * x[N - 1 - j];
				result[j] += w * U[M * i + j];
			}
			
		}
	}
	else if (option == 'l')
	{
                #pragma omp parallel for
		for (int i = 0; i < M; i++)
			result[i] = 0.0;
		for (int i = 0; i < R; i++)
		{
			w = 0.0;
			for (int j = 0; j < M; j++)
			{
				if (j < N)
					w += V[N * i + j] * x[j];
				result[j] += w * U[M * i + j];
			}
		}
	}
	else if (option == 'u')
	{
                #pragma omp parallel for
		for (int i = 0; i < M; i++)
			result[i] = 0.0;
		for (int i = 0; i < R; i++)
		{
			w = 0.0;
			for (int j = 0; j < M; j++)
			{
				if (j < N)
					w += V[N * i + (N - 1 - j)] * x[N - 1 - j];

				result[N - 1 - j] += w * U[M * i + N - 1 - j];
			}
			
		}
	}
	else if (option == 'q')
	{
                #pragma omp parallel for
		for (int i = 0; i < M; i++)
			result[i] = 0.0;
		for (int i = 0; i < R; i++)
		{
			w = 0.0;
			for (int j = 0; j < M; j++)
			{
				if (j < N)
					w += V[N * i + j] * x[j];

				result[N - j - 1] += w * U[M * i + N - j - 1];
			}
		}

	}
	else
	{
		printf("incorrect matvec option\nreturn null\n");
		free(result);
		result = NULL;
	}
	return result;
}
#endif
void TCross_Parallel_v1::get_approximations(double *outputU, double *outputV)
{
    int R = this->get_rank();
    int N = this->rows_number;
    int M = this->columns_number;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < R; j++) *(outputU+j+R*i) = U[j+R*i];
    for (int i = 0; i < M; i++)
        for (int j = 0; j < R; j++) *(outputV+j+R*i) = U[j+R*i];

}

TCross_Parallel_v1::~TCross_Parallel_v1()
{
    //*******************************RISHAT*ADDED*BEGIN****************************************
#ifdef CUDA_FFT
    dealloc_cuda_data();
#endif
    //*******************************RISHAT*ADDED*END******************************************
    free(U);
    free(V);
    free(C);
    free(hat_A_inv);
    free(RT);
    free(CAT);
    free(AR);
    free(rows_numbers);
    free(columns_numbers);
}
