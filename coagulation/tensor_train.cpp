#include "tensor_train.h"
#include <iostream>
using namespace std;

MulTTensor::MulTTensor():TTensor()
{
}

MulTTensor::MulTTensor(const int &d, int *&m):TTensor(d, m)
{
}

void MulTTensor::set_trains( const TTensorTrain *first,  const TTensorTrain *second)
{
	tt_one = const_cast <TTensorTrain*>(first);
	tt_two = const_cast <TTensorTrain*>(second);
}

double MulTTensor::operator[](int *i)
{
	return (this->tt_one[0][i]) * (this->tt_two[0][i]) ; 
}

TTensorTrain::TTensorTrain(const TTensorTrain &tens, const  double &c)
{
    this->dimensionality = tens.dimensionality;
    carriages = (double **) malloc(this->dimensionality*sizeof(double*));
    this->modes_sizes = (int *) malloc(this->dimensionality * sizeof(int));
    for (int i = 0; i < this->dimensionality; i++)
    {
        this->modes_sizes[i] = tens.modes_sizes[i];
    }
    this->ranks = (int *) malloc((this->dimensionality + 1) * sizeof(int));
    for (int i = 0; i <= this->dimensionality; i++)
    {
        this->ranks[i] = 1;
    }
    for (int k = 0; k < this->dimensionality; k++)
    {
        this->carriages[k] = (double *) malloc(this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1] * sizeof(double));
        for (int i = 0; i < this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1]; i++)
        {
            this->carriages[k][i] = c;
        }
    }
}

TTensorTrain::TTensorTrain(const TTensorTrain &tens, const  double *c)
{
    this->dimensionality = tens.dimensionality;
    carriages = (double **) malloc(this->dimensionality*sizeof(double*));
    this->modes_sizes = (int *) malloc(this->dimensionality * sizeof(int));
    for (int i = 0; i < this->dimensionality; i++)
    {
        this->modes_sizes[i] = tens.modes_sizes[i];
    }
    this->ranks = (int *) malloc((this->dimensionality + 1) * sizeof(int));
    for (int i = 0; i <= this->dimensionality; i++)
    {
        this->ranks[i] = 1;
    }
    for (int k = 0; k < this->dimensionality; k++)
    {
        this->carriages[k] = (double *) malloc(this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1] * sizeof(double));
        for (int i = 0; i < this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1]; i++)
        {
            this->carriages[k][i] = c[i];
        }
    }
}

        
/*
void TTensorTrain::Orthogonalize()
{
    for (int k = 0; k < this->dimensionality - 1; k++)
    {
        double *U, *VT;
        orthofactors(ranks[k + 1], ranks[k] * this->modes_sizes[k], carriages[k], U, VT);
        free(carriages[k]);
        carriages[k] = VT;
        for (int i = 0; i < this->modes_sizes[k + 1]; i++)
        {
            double tmp[ranks[k + 2] * ranks[k + 1]];
            gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ranks[k + 2], ranks[k + 1], ranks[k + 1], 1e0, carriages[k + 1] + i * ranks[k + 2] * ranks[k + 1], ranks[k + 2], U, ranks[k + 1], 0e0, tmp, ranks[k + 2]);
            copy(ranks[k + 2] * ranks[k + 1], tmp, 1, carriages[k + 1] + i * ranks[k + 2] * ranks[k + 1], 1);
        }
        free(U);
    }
}


void TTensorTrain::Compress(const double &tol)
{
    for (int k = 0; k < this->dimensionality - 1; k++)
    {
        double *U = (double *) malloc(ranks[k + 1] * ranks[k] * this->modes_sizes[k] * sizeof(double));
        double *V = (double *) malloc(ranks[k + 2] * ranks[k + 1] * this->modes_sizes[k + 1] * sizeof(double));
        double alpha[ranks[k + 1]], beta[ranks[k + 1] - 1], W[ranks[k + 1] * ranks[k + 1]], buf[ranks[k + 1]];
        if (true)//(ranks[k] * this->modes_sizes[k] <= ranks[k + 2] * this->modes_sizes[k + 1])
        {
            for (int i = 0; i < ranks[k] * this->modes_sizes[k]; i++)
            {
                random(U[i]);
            }
            scal(ranks[k] * this->modes_sizes[k], ((double) 1e0)/nrm2(ranks[k] * this->modes_sizes[k], U, 1), U, 1);
            double tmp[ranks[k] * this->modes_sizes[k]];
            for (int i = 0; i < ranks[k] * this->modes_sizes[k]; i++)
            {
                tmp[i] = conj(U[i]);
            }
            gemv(CblasColMajor, CblasNoTrans, ranks[k + 1], ranks[k] * this->modes_sizes[k], 1e0, carriages[k], ranks[k + 1], tmp, 1, 0e0, buf, 1);
            for (int i = 0; i < this->modes_sizes[k + 1]; i++)
            {
                gemv(CblasColMajor, CblasNoTrans, ranks[k + 2], ranks[k + 1], 1e0, carriages[k + 1] + i * ranks[k + 2] * ranks[k + 1], ranks[k + 2], buf, 1, 0e0, V + i * ranks[k + 2], 1);
            }
            alpha[0] = nrm2(ranks[k + 2] * this->modes_sizes[k + 1], V, 1);
            scal(ranks[k + 2] * this->modes_sizes[k + 1], ((double) 1e0)/alpha[0], V, 1);
            W[0] = (double) 1e0;
            for (int i = 1; i < ranks[k + 1]; i++)
            {
                copy(ranks[k] * this->modes_sizes[k], U + (i-1) * ranks[k] * this->modes_sizes[k], 1, U + i * ranks[k] * this->modes_sizes[k], 1);
                for (int j = 0; j < ranks[k + 1]; j++) buf[j] = 0e0;
                for (int j = 0; j < this->modes_sizes[k + 1]; j++)
                {
                    gemv(CblasTrans, ranks[k + 1], ranks[k + 2], 1e0, carriages[k + 1] + j * ranks[k + 2] * ranks[k + 1], ranks[k + 2], V + (i - 1) * ranks[k + 2] * (this->modes_sizes[k + 1] + j), 1, 1e0, buf, 1);
                }
                gemv(CblasTrans, ranks[k] * this->modes_sizes[k], ranks[k + 1], 1e0, carriages[k], ranks[k + 1], buf, 1, - alpha[i - 1], U + i * ranks[k] * this->modes_sizes[k], 1);
                beta[i - 1] = nrm2(ranks[k] * this->modes_sizes[k], U + i * ranks[k] * this->modes_sizes[k], 1);
                scal(ranks[k] * this->modes_sizes[k], 1e0/beta[i - 1], U + i * ranks[k] * this->modes_sizes[k], 1);
                copy(ranks[k + 2] * this->modes_sizes[k + 1], V + (i - 1) * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                for (int j = 0; j < ranks[k] * this->modes_sizes[k]; j++)
                {
                    tmp[j] = conj(U[j + i * ranks[k] * this->modes_sizes[k]]);
                }
                gemv(CblasColMajor, CblasNoTrans, ranks[k + 1], ranks[k] * this->modes_sizes[k], 1e0, carriages[k], ranks[k + 1], tmp, 1, 0e0, buf, 1);
                for (int j = 0; j < this->modes_sizes[k + 1]; j++)
                {
                    gemv(CblasColMajor, CblasNoTrans, ranks[k + 2], ranks[k + 1], 1e0, carriages[k + 1] + j * ranks[k + 2] * ranks[k + 1], ranks[k + 2], buf, 1, 0e0, V + ranks[k + 2] * (i * this->modes_sizes[k + 1] + j), 1);
                }
                alpha[i] = nrm2(ranks[k + 2] * this->modes_sizes[k + 1], V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                scal(ranks[k + 2] * this->modes_sizes[k + 1], ((double) 1e0)/alpha[i], V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                W[i + i * ranks[k + 1]] = 1e0;
                for (int j = 0; j < i; j++)
                {
                    W(i + j * ranks[k + 1]) = dot(ranks[k] * this->modes_sizes[k], U + i * ranks[k] * this->modes_sizes[k], 1, U + j * ranks[k] * this->modes_sizes[k], 1);
                    W(j + i * ranks[k + 1]) = dot(ranks[k + 2] * this->modes_sizes[k + 1], V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + j * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                }
                if (abs(W(i + ranks[k + 1] * iamax(i, W + i, ranks[k + 1]))) > sqrt(eps))
                {
                    for (int j = 0; j < i - 1; j++)
                    {
                        axpy(ranks[k] * this->modes_sizes[k], - W[i-1 + j * ranks[k + 1]], U + j * ranks[k] * this->modes_sizes[k], 1, U + (i-1) * ranks[k] * this->modes_sizes[k], 1);
                        axpy(ranks[k] * this->modes_sizes[k], - W[i + j * ranks[k + 1]], U + j * ranks[k] * this->modes_sizes[k], 1, U + i * ranks[k] * this->modes_sizes[k], 1);
                        W[i  + j * ranks[k + 1]] = dot(ranks[k] * this->modes_sizes[k], U + i * ranks[k] * this->modes_sizes[k], 1, U + j * ranks[k] * this->modes_sizes[k], 1);
                    }
                    scal(ranks[k] * this->modes_sizes[k], 1e0/nrm2(ranks[k] * this->modes_sizes[k], U + (i-1) * ranks[k] * this->modes_sizes[k], 1), U + (i-1) * ranks[k] * this->modes_sizes[k], 1);
                    axpy(ranks[k] * this->modes_sizes[k], - W[i + (i - 1) * ranks[k + 1]], U + (i - 1) * ranks[k] * this->modes_sizes[k], 1, U + i * ranks[k] * this->modes_sizes[k], 1);
                    scal(ranks[k] * this->modes_sizes[k], 1e0/nrm2(ranks[k] * this->modes_sizes[k], U + i * ranks[k] * this->modes_sizes[k], 1), U + i * ranks[k] * this->modes_sizes[k], 1);
                    W[i  + (i - 1) * ranks[k + 1]] = dot(ranks[k] * this->modes_sizes[k], U + i * ranks[k] * this->modes_sizes[k], 1, U + (i - 1) * ranks[k] * this->modes_sizes[k], 1);
                }
                if (abs(W(i * ranks[k + 1] + iamax(i, W + i * ranks[k + 1], 1))) > sqrt(eps))
                {
                    for (int j = 0; j < i - 1; j++)
                    {
                        axpy(ranks[k + 2] * this->modes_sizes[k + 1], - W[(i-1) * ranks[k + 1] + j], V + j * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + (i-1) * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                        axpy(ranks[k + 2] * this->modes_sizes[k + 1], - W[i * ranks[k + 1] + j], V + j * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                        W[i * ranks[k + 1]  + j) = dot(ranks[k + 2] * this->modes_sizes[k + 1], V + j * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + j * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                    }
                    scal(ranks[k + 2] * this->modes_sizes[k + 1], 1e0/nrm2(ranks[k] * this->modes_sizes[k],  V + (i-1) * ranks[k + 2] * this->modes_sizes[k + 1], 1), V + (i-1) * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                    axpy(ranks[k + 2] * this->modes_sizes[k + 1], - W[i * ranks[k + 1] + i - 1], V + (i - 1) * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                    scal(ranks[k + 2] * this->modes_sizes[k + 1], 1e0/nrm2(ranks[k + 2] * this->modes_sizes[k + 1], V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1), V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                    W[i * ranks[k + 1]  + i - 1) = dot(ranks[k + 2] * this->modes_sizes[k + 1], V + i * ranks[k + 2] * this->modes_sizes[k + 1], 1, V + (i - 1) * ranks[k + 2] * this->modes_sizes[k + 1], 1);
                }
            }
            free(U);
            free(V);
        }
    }
    for (int k = 0; k < this->dimensionality - 1; k++)
    {
        double *U, *VT;
        int newrank;
        reducedorthofactors(tol / sqrt((double) this->dimensionality - 1), ranks[k + 1], ranks[k] * this->modes_sizes[k], carriages[k], U, VT, newrank);
        free(carriages[k]);
        carriages[k] = VT;
        double *tmp = (double *) malloc(ranks[k + 2] * newrank * this->modes_sizes[k + 1] * sizeof(double));
        for (int i = 0; i < this->modes_sizes[k + 1]; i++)
        {
            gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ranks[k + 2], newrank, ranks[k + 1], 1e0, carriages[k + 1] + i * ranks[k + 2] * ranks[k + 1], ranks[k + 2], U, ranks[k + 1], 0e0, tmp + i * ranks[k + 2] * newrank, ranks[k + 2]);
        }
        ranks[k + 1] = newrank;
        free(carriages[k + 1]);
        carriages[k + 1] = tmp;
        free(U);
    }
}*/

void TTensorTrain::save_to_file(char *filename)
{
    FILE *f;
    f = fopen(filename, "wb");
    fwrite(&this->dimensionality, sizeof(int), 1, f);
    fwrite(this->modes_sizes, sizeof(int), this->dimensionality, f);
    fwrite(this->ranks, sizeof(int), this->dimensionality + 1, f);
    for (int k = 0; k < this->dimensionality; k++)
    {
        fwrite(this->carriages[k], sizeof(double), this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1], f);
    }
    fclose(f);
}

void TTensorTrain::load_from_file(char *filename)
{
    FILE *f;
    f = fopen(filename, "rb");
    for (int k = 0; k < this->dimensionality; k++)
    {
        free(carriages[k]);
    }
    free(carriages);
    carriages = NULL;
    free(ranks);
    ranks = NULL;
    fread(&this->dimensionality, sizeof(int), 1, f);
    this->modes_sizes = (int *) realloc(this->modes_sizes, this->dimensionality*sizeof(int));
    ranks = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    carriages = (double **) malloc(this->dimensionality*sizeof(double*));
    fread(this->modes_sizes, sizeof(int), this->dimensionality, f);
    fread(this->ranks, sizeof(int), this->dimensionality + 1, f);
    for (int k = 0; k < this->dimensionality; k++)
    {
        carriages[k] = (double *) malloc(this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1] * sizeof(double));
        fread(this->carriages[k], sizeof(double), this->ranks[k] * this->modes_sizes[k] * this->ranks[k + 1], f);
    }
    fclose(f);
}

double tt_dot(int d, int * modes_sizes, int * a_ranks, double ** a, int * b_ranks, double ** b)
{
    double *v = (double *) malloc(sizeof(double));
    v[0] = 1.0;
    for (int k = 0; k < d; k++)
    {
        double *g = (double *) malloc(a_ranks[k]*b_ranks[k + 1]*sizeof(double));
        double *v1 = (double *) calloc(a_ranks[k + 1]*b_ranks[k + 1], sizeof(double));
        for (int i = 0; i < modes_sizes[k]; i++)
        {
            gemm(CblasColMajor, CblasNoTrans, CblasConjTrans, a_ranks[k], b_ranks[k + 1], b_ranks[k], 1.0, v, a_ranks[k], b[k] + i*b_ranks[k]*b_ranks[k + 1], b_ranks[k + 1], 0.0, g, a_ranks[k]);
            gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a_ranks[k + 1], b_ranks[k + 1], a_ranks[k], 1.0, a[k] + i*a_ranks[k]*a_ranks[k + 1], a_ranks[k + 1], g, a_ranks[k], 1.0, v1, a_ranks[k + 1]);
        }
        free(v);
        free(g);
        v = v1;
    }
    double r = v[0];
    free(v);
    return r;
}

double TTensorTrain::nrm2()
{
    return sqrt(abs(tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, ranks, carriages)));
}

double TTensorTrain::operator[](int *i)
{
    double *v = (double *) malloc(sizeof(double));
    v[0] = 1.0;
    for (int k = 0; k < this->dimensionality; k++)
    {
        double *v1 = (double *) malloc(ranks[k + 1] * sizeof(double));
        gemv(CblasColMajor, CblasNoTrans, ranks[k + 1], ranks[k], 1.0, carriages[k] + i[k] * ranks[k] * ranks[k + 1], ranks[k + 1], v, 1, 0.0, v1, 1);
        free(v);
        v = v1;
    }
    double result = v[0];
    free(v);
    return result;
}

double TTensorTrain :: operator[](const int &ind) 
{
    int * multindex = (int *) malloc((dimensionality - 1) * sizeof(int));
    int index = ind;
    for (int i = 0; i < dimensionality; i++){
        multindex[i] = index % modes_sizes[i];
        index /= modes_sizes[i];
        // it's correct
    }
    return (*this)[multindex];
}
TTensorTrain::~TTensorTrain()
{
    for (int k = 0; k < this->get_dimensionality(); k++)
    {
        free(carriages[k]);
    }
    free(carriages);
    carriages = NULL;
    free(ranks);
    ranks = NULL;
}

TTensorTrain::TTensorTrain()
{
    ranks = NULL;
    carriages = NULL;
}

void TTensorTrain::Approximate(TTensor * tensor, const TTensorTrainParameters & parameters)
{
    for (int k = 0; k < this->dimensionality; k++)
    {
        free(carriages[k]);
    }
    free(carriages);
    carriages = NULL;
    free(ranks);
    ranks = NULL;
    this->dimensionality = tensor->get_dimensionality();
    this->modes_sizes = (int *) realloc(this->modes_sizes, this->dimensionality*sizeof(int));
    int max_mode_size = 0;
    for (int i = 0; i < this->dimensionality; i++)
    {
        this->modes_sizes[i] = tensor->get_mode_size(i);
        if (this->modes_sizes[i] > max_mode_size) max_mode_size = this->modes_sizes[i];
    }
    int ***left_indeces = (int ***) malloc((this->dimensionality + 1)*sizeof(int**));
    int ***right_indeces = (int ***) malloc((this->dimensionality + 1)*sizeof(int**));
    int *left_indeces_number = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    int *right_indeces_number = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    ranks = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    carriages = (double **) malloc(this->dimensionality*sizeof(double*));
    for (int i = 0; i <= this->dimensionality; i++) ranks[i] = 0;
    left_indeces_number[0] = 0;
    left_indeces_number[this->dimensionality] = 0;
    left_indeces[0] = (int**) calloc(left_indeces_number[0], sizeof(int*));
    left_indeces[this->dimensionality] = NULL;
    right_indeces_number[0] = 0;
    right_indeces_number[this->dimensionality] = 0;
    right_indeces[this->dimensionality] = (int**) calloc(right_indeces_number[this->dimensionality], sizeof(int*));
    right_indeces[0] = NULL;
    for (int i = 1; i < this->dimensionality; i++)
    {
        left_indeces_number[i] = 0;
        left_indeces[i] = NULL;
        right_indeces_number[i] = 0;
        right_indeces[i] = NULL;
    }
    int ***back_left_indeces = (int ***) malloc((this->dimensionality + 1)*sizeof(int**));
    int ***back_right_indeces = (int ***) malloc((this->dimensionality + 1)*sizeof(int**));
    int *back_left_indeces_number = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    int *back_right_indeces_number = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    int *back_ranks = (int *) malloc((this->dimensionality + 1)*sizeof(int));
    double **back_carriages = (double **) malloc(this->dimensionality*sizeof(double*));
    for (int i = 0; i <= this->dimensionality; i++) back_ranks[i] = 0;
    back_left_indeces_number[0] = 0;
    back_left_indeces_number[this->dimensionality] = 0;
    back_left_indeces[0] = (int**) calloc(back_left_indeces_number[0], sizeof(int*));
    back_left_indeces[this->dimensionality] = NULL;
    back_right_indeces_number[0] = 0;
    back_right_indeces_number[this->dimensionality] = 0;
    back_right_indeces[this->dimensionality] = (int**) calloc(back_right_indeces_number[this->dimensionality], sizeof(int*));
    back_right_indeces[0] = NULL;
    for (int i = 1; i < this->dimensionality; i++)
    {
        back_left_indeces_number[i] = 0;
        back_left_indeces[i] = NULL;
        back_right_indeces_number[i] = 0;
        back_right_indeces[i] = NULL;
    }
    for (int i = 0; i < this->dimensionality; i++)
    {
        carriages[i] = NULL;
        back_carriages[i] = NULL;
    }
    TCross_Parallel_v1 cross_approx;
    TCross_Parallel_v1_Parameters param;
    param.tolerance = (parameters.tolerance / sqrt((double) this->dimensionality));
    param.maximal_iterations_number = 0;//parameters.cross_max_iterations;
    param.number_of_checked_elements = 1;
    param.stop_rank = 0;//parameters.stop_rank;
    param.max_rank =10;// max(parameters.stop_rank, 1);
    int iters = 0;
    while ((iters < parameters.maximal_iterations_number) || (parameters.maximal_iterations_number == 0))
    {
        //cout << "Sweep " << iters << endl;
        ranks[0] = 1;
        ranks[this->dimensionality] = 1;
        param.start_from_column = true;
        for (int k = 1; k < this->dimensionality; k++)
        {
            //cout << k;
            for (int o = 0; o < left_indeces_number[k]; o++) free(left_indeces[k][o]);
            free(left_indeces[k]);
            left_indeces_number[k] = ranks[k - 1] * this->get_mode_size(k - 1);
//            left_indeces[k] = (int **) malloc(left_indeces_number[k]*sizeof(int*));
            left_indeces[k] = new int*[left_indeces_number[k]];
            for (int i = 0; i < this->get_mode_size(k - 1); i++)
            {
                for (int j = 0; j < ranks[k - 1]; j++)
                {
                    left_indeces[k][i*ranks[k - 1]+j] = (int *) malloc(k*sizeof(int));
                    for (int o = 0; o < k - 1; o++)
                    {
                        left_indeces[k][i*ranks[k - 1]+j][o] = left_indeces[k - 1][j][o];
                    }
                    left_indeces[k][i*ranks[k - 1]+j][k - 1] = i;
                }
            }
            if (k != this->dimensionality - 1)
            {
                right_indeces[k] = (int **) realloc(right_indeces[k], max(right_indeces_number[k] + back_right_indeces_number[k] + right_indeces_number[k - 1], left_indeces_number[k])*sizeof(int *));
                for (int i = 0; i < back_right_indeces_number[k]; i++)
                {
                    right_indeces[k][right_indeces_number[k] + i] = (int *) malloc((this->dimensionality - k)*sizeof(int));
                    for (int o = 0; o < this->dimensionality - k; o++)
                    {
                        right_indeces[k][right_indeces_number[k] + i][o] = back_right_indeces[k][i][o];
                    }
                }
                right_indeces_number[k] += back_right_indeces_number[k];
                for (int i = 0; i < right_indeces_number[k - 1]; i++)
                {
                    right_indeces[k][right_indeces_number[k] + i] = (int *) malloc((this->dimensionality - k)*sizeof(int));
                    for (int o = 0; o < this->dimensionality - k; o++)
                    {
                        right_indeces[k][right_indeces_number[k] + i][o] = right_indeces[k - 1][i][o];
                    }
                }
                right_indeces_number[k] += right_indeces_number[k - 1];
                if (right_indeces_number[k] < left_indeces_number[k])
                {
                    for (int i = right_indeces_number[k]; i < left_indeces_number[k]; i++)
                    {
                        right_indeces[k][i] = (int *) malloc((this->dimensionality - k)*sizeof(int));
                        for (int o = 0; o < this->dimensionality - k; o++)
                            right_indeces[k][i][o] = rand()%(this->get_mode_size(this->dimensionality - o - 1));
                    }
                    right_indeces_number[k] = left_indeces_number[k];
                }
            }
            else
            {
                for (int i = 0; i < right_indeces_number[k]; i++)
                {
                    free(right_indeces[k][i]);
                }
                right_indeces_number[k] = this->get_mode_size(this->dimensionality - 1);
                right_indeces[k] = (int **) realloc(right_indeces[k], right_indeces_number[k] * sizeof(int*));
                for (int i = 0; i < right_indeces_number[k];  i++)
                {
                    right_indeces[k][i] = (int *) malloc(sizeof(int));
                    right_indeces[k][i][0] = i;
                }
            }
            TUnfoldingSubMatrix work_matrix(tensor, k, left_indeces_number[k], right_indeces_number[k], left_indeces[k], right_indeces[k]);
            param.max_rank = ranks[k];
            cross_approx.Approximate(&work_matrix, param);
            ranks[k] = cross_approx.get_rank();
            //cout << " rank = " << ranks[k] << endl;
            carriages[k - 1] = (double *) malloc(ranks[k - 1] * ranks[k] * this->get_mode_size(k - 1) * sizeof(double));
            copy(ranks[k - 1] * ranks[k] * this->get_mode_size(k - 1), cross_approx.export_CAT(), 1, carriages[k - 1], 1);
//            gemm(CblasColMajor, CblasTrans, CblasTrans, ranks[k], ranks[k - 1] * this->get_mode_size(k - 1), ranks[k], 1.0, cross_approx.export_hat_A_inv(), ranks[k], cross_approx.export_C(), ranks[k - 1] * this->get_mode_size(k - 1), 0.0, carriages[k - 1], ranks[k]);
            int ** l_indeces = (int **) malloc(ranks[k]*sizeof(int *)), ** r_indeces = (int **) malloc(ranks[k]*sizeof(int *));
            for (int i = 0; i < ranks[k]; i++)
            {
                l_indeces[i] = left_indeces[k][cross_approx.get_row_number(i)];
                left_indeces[k][cross_approx.get_row_number(i)] = NULL;
                r_indeces[i] = right_indeces[k][cross_approx.get_column_number(i)];
                right_indeces[k][cross_approx.get_column_number(i)] = NULL;
            }
            for (int o = 0; o < left_indeces_number[k]; o++) free(left_indeces[k][o]);
            for (int o = 0; o < right_indeces_number[k]; o++) free(right_indeces[k][o]);
            free(left_indeces[k]);
            free(right_indeces[k]);
            left_indeces_number[k] = ranks[k];
            right_indeces_number[k] = ranks[k];
            left_indeces[k] = l_indeces;
            right_indeces[k] = r_indeces;
            l_indeces = NULL;
            r_indeces = NULL;
        }
        carriages[this->dimensionality - 1] = (double *) malloc(ranks[this->dimensionality - 1]*ranks[this->dimensionality]*this->get_mode_size(this->dimensionality - 1)*sizeof(double));
        for (int i = 0; i < ranks[this->dimensionality - 1]*ranks[this->dimensionality]; i++)
        {
            copy(this->get_mode_size(this->dimensionality - 1), cross_approx.export_RT() + i*this->get_mode_size(this->dimensionality - 1), 1, carriages[this->dimensionality - 1] + i, ranks[this->dimensionality - 1]*ranks[this->dimensionality]);
        }
        if (iters != 0)
        {
            double last_norm = tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, ranks, carriages), diff_norm = last_norm + tt_dot(this->dimensionality, this->modes_sizes, back_ranks, back_carriages, back_ranks, back_carriages) - 2.0 * tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, back_ranks, back_carriages);
            //cout << "Diff = " << sqrt(diff_norm / last_norm) << endl;
            if(real(diff_norm) <= real(parameters.tolerance * parameters.tolerance * last_norm))
            {
                for (int k = 0; k < this->dimensionality; k++)
                {
                    free(back_carriages[k]);
                    back_carriages[k] = NULL;
                }
                for (int i = 0; i <= this->dimensionality; i++)
                {
                    for (int o = 0; o < left_indeces_number[i]; o++) free(left_indeces[i][o]);
                    for (int o = 0; o < right_indeces_number[i]; o++) free(right_indeces[i][o]);
                    for (int o = 0; o < back_left_indeces_number[i]; o++) free(back_left_indeces[i][o]);
                    for (int o = 0; o < back_right_indeces_number[i]; o++) free(back_right_indeces[i][o]);
                    free(right_indeces[i]);
                    free(left_indeces[i]);
                    free(back_left_indeces[i]);
                    free(back_right_indeces[i]);
                }
                free(right_indeces);
                free(left_indeces);
                free(back_left_indeces);
                free(back_right_indeces);
                free(right_indeces_number);
                free(left_indeces_number);
                free(back_left_indeces_number);
                free(back_right_indeces_number);
                free(back_ranks);
                free(back_carriages);
                return;
            }
        }
        //cout << "Back sweep " << iters << endl;;
        for (int k = 0; k < this->dimensionality; k++)
        {
            free(back_carriages[k]);
            back_carriages[k] = NULL;
        }
        back_ranks[0] = 1;
        back_ranks[this->dimensionality] = 1;
        param.start_from_column = false;
        for (int k = this->dimensionality - 1; k > 0; k--)
        {
            //cout << k;
            for (int o = 0; o < back_right_indeces_number[k]; o++) free(back_right_indeces[k][o]);
            free(back_right_indeces[k]);
            back_right_indeces_number[k] = back_ranks[k + 1] * this->get_mode_size(k);
            back_right_indeces[k] = (int **) malloc(back_right_indeces_number[k]*sizeof(int *));
            for (int i = 0; i < this->get_mode_size(k); i++)
            {
                for (int j = 0; j < back_ranks[k + 1]; j++)
                {
                    back_right_indeces[k][i*back_ranks[k + 1]+j] = (int *) malloc((this->dimensionality - k) * sizeof(int));
                    for (int o = 0; o < this->dimensionality - k - 1; o++)
                        back_right_indeces[k][i*back_ranks[k + 1]+j][o] = back_right_indeces[k + 1][j][o];
                    back_right_indeces[k][i*back_ranks[k + 1]+j][this->dimensionality - k - 1] = i;
                }
            }
            if (k != 1)
            {
                back_left_indeces[k] = (int **) realloc(back_left_indeces[k], max(back_left_indeces_number[k] + left_indeces_number[k] + back_left_indeces_number[k + 1], back_right_indeces_number[k])*sizeof(int*));
                for (int i = 0; i < left_indeces_number[k]; i++)
                {
                    back_left_indeces[k][back_left_indeces_number[k] + i] = (int *) malloc(k * sizeof(int));
                    for (int o = 0; o < k; o++)
                        back_left_indeces[k][back_left_indeces_number[k] + i][o] = left_indeces[k][i][o];
                }
                back_left_indeces_number[k] += left_indeces_number[k];
                for (int i = 0; i < back_left_indeces_number[k + 1]; i++)
                {
                    back_left_indeces[k][back_left_indeces_number[k] + i] = (int *) malloc(k * sizeof(int));
                    for (int o = 0; o < k; o++)
                        back_left_indeces[k][back_left_indeces_number[k] + i][o] = back_left_indeces[k + 1][i][o];
                }
                back_left_indeces_number[k] += back_left_indeces_number[k + 1];
                if (back_left_indeces_number[k] < back_right_indeces_number[k])
                {
                    for (int i = back_left_indeces_number[k]; i < back_right_indeces_number[k]; i++)
                    {
                        back_left_indeces[k][i] = (int *) malloc(k * sizeof(int));
                        for (int o = 0; o < k; o++)
                            back_left_indeces[k][i][o] = rand()%this->get_mode_size(o);
                    }
                    back_left_indeces_number[k] = back_right_indeces_number[k];
                }
            }
            else
            {
                for (int i = 0; i < back_left_indeces_number[k]; i++)
                {
                    free(back_left_indeces[k][i]);
                }
                back_left_indeces_number[k] = this->get_mode_size(0);
                back_left_indeces[k] = (int **) realloc(back_left_indeces[k], back_left_indeces_number[k] * sizeof(int*));
                for (int i = 0; i < back_left_indeces_number[k]; i++)
                {
                    back_left_indeces[k][i] = (int *) malloc(sizeof(int));
                    back_left_indeces[k][i][0] = i;
                }
            }
            TUnfoldingSubMatrix work_matrix(tensor, k, back_left_indeces_number[k], back_right_indeces_number[k], back_left_indeces[k], back_right_indeces[k]);
            cross_approx.Approximate(&work_matrix, param);
            back_ranks[k] = cross_approx.get_rank();
            //cout << " rank = " << back_ranks[k] << endl;
            back_carriages[k] = (double *) malloc(back_ranks[k + 1] * back_ranks[k] * this->get_mode_size(k) * sizeof(double));
            for (int i = 0; i < this->get_mode_size(k); i++)
            {
//                gemm(CblasColMajor, CblasNoTrans, CblasTrans, back_ranks[k + 1], back_ranks[k], back_ranks[k], 1.0, cross_approx.export_RT() + i*back_ranks[k + 1], back_ranks[k + 1] * this->get_mode_size(k), cross_approx.export_hat_A_inv(), back_ranks[k], 0.0, back_carriages[k] + i * back_ranks[k] * back_ranks[k + 1], back_ranks[k + 1]);
              for(int j = 0; j < back_ranks[k + 1]; j++)
                {
                    copy(back_ranks[k], cross_approx.export_AR() + back_ranks[k] * (j + back_ranks[k + 1] * i), 1, back_carriages[k] + j + back_ranks[k] * back_ranks[k + 1] * i, back_ranks[k + 1]);
                }
            }
            
            int ** l_indeces = (int **) malloc(back_ranks[k]*sizeof(int*)), ** r_indeces = (int **) malloc(back_ranks[k]*sizeof(int*));
            for (int i = 0; i < back_ranks[k]; i++)
            {
                l_indeces[i] = back_left_indeces[k][cross_approx.get_row_number(i)];
                back_left_indeces[k][cross_approx.get_row_number(i)] = NULL;
                r_indeces[i] = back_right_indeces[k][cross_approx.get_column_number(i)];
                back_right_indeces[k][cross_approx.get_column_number(i)] = NULL;
            }
            for (int o = 0; o < back_left_indeces_number[k]; o++) free(back_left_indeces[k][o]);
            for (int o = 0; o < back_right_indeces_number[k]; o++) free(back_right_indeces[k][o]);
            free(back_left_indeces[k]);
            free(back_right_indeces[k]);
            back_left_indeces_number[k] = back_ranks[k];
            back_right_indeces_number[k] = back_ranks[k];
            back_left_indeces[k] = l_indeces;
            back_right_indeces[k] = r_indeces;
            l_indeces = NULL;
            r_indeces = NULL;
        }
        back_carriages[0] = (double *) malloc(back_ranks[1] * back_ranks[0] * this->get_mode_size(0) * sizeof(double));
        for (int i = 0; i < back_ranks[0]*back_ranks[1]; i++)
        {
            copy(this->get_mode_size(0), cross_approx.export_C() + i*this->get_mode_size(0), 1, back_carriages[0] + i, back_ranks[0]*back_ranks[1]);
        }
        double last_norm = tt_dot(this->dimensionality, this->modes_sizes, back_ranks, back_carriages, back_ranks, back_carriages), diff_norm = last_norm + tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, ranks, carriages) - 2.0 * tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, back_ranks, back_carriages);
        //cout << "Diff = " << sqrt(diff_norm / last_norm) << endl;
        if(abs(real(diff_norm)) <= abs(real(parameters.tolerance * parameters.tolerance * last_norm)))
        {
            for (int k = 0; k < this->dimensionality; k++)
            {
                free(carriages[k]);
                carriages[k] = NULL;
            }
            for (int i = 0; i <= this->dimensionality; i++)
            {
                for (int o = 0; o < left_indeces_number[i]; o++) free(left_indeces[i][o]);
                for (int o = 0; o < right_indeces_number[i]; o++) free(right_indeces[i][o]);
                for (int o = 0; o < back_left_indeces_number[i]; o++) free(back_left_indeces[i][o]);
                for (int o = 0; o < back_right_indeces_number[i]; o++) free(back_right_indeces[i][o]);
                free(right_indeces[i]);
                free(left_indeces[i]);
                free(back_left_indeces[i]);
                free(back_right_indeces[i]);
            }
            free(right_indeces);
            free(left_indeces);
            free(back_left_indeces);
            free(back_right_indeces);
            free(right_indeces_number);
            free(left_indeces_number);
            free(back_left_indeces_number);
            free(back_right_indeces_number);
            free(ranks);
            free(carriages);
            ranks = back_ranks;
            carriages = back_carriages;
            return;
        }
        for (int k = 0; k < this->dimensionality; k++)
        {
            free(carriages[k]);
            carriages[k] = NULL;
        }
        iters++;
    }
    for (int i = 0; i <= this->dimensionality; i++)
    {
        for (int o = 0; o < left_indeces_number[i]; o++) free(left_indeces[i][o]);
        for (int o = 0; o < right_indeces_number[i]; o++) free(right_indeces[i][o]);
        for (int o = 0; o < back_left_indeces_number[i]; o++) free(back_left_indeces[i][o]);
        for (int o = 0; o < back_right_indeces_number[i]; o++) free(back_right_indeces[i][o]);
        free(right_indeces[i]);
        free(left_indeces[i]);
        free(back_left_indeces[i]);
        free(back_right_indeces[i]);
    }
    free(right_indeces);
    free(left_indeces);
    free(back_left_indeces);
    free(back_right_indeces);
    free(right_indeces_number);
    free(left_indeces_number);
    free(back_left_indeces_number);
    free(back_right_indeces_number);
    free(ranks);
    free(carriages);
    ranks = back_ranks;
    carriages = back_carriages;
}


int TTensorTrain::get_rank(int &i)
{
    return ranks[i];
}

double TTensorTrain::dot(const TTensorTrain &tens)
{
    return tt_dot(this->dimensionality, this->modes_sizes, ranks, carriages, tens.ranks, tens.carriages);
}

void TTensorTrain :: SVDCompress (const double eps, int cutrank)
{

	int maxrank = this->ranks[0];
	int maxcarsize = 0;

	for (int i = 0; i < this->dimensionality ; i++) 
	{
		maxrank = max(maxrank, this->ranks[i + 1]);
		maxcarsize = max(maxcarsize, this->ranks[i + 1] * this->ranks[i] * this->modes_sizes[i]);
	}

	if (cutrank == 0) 
	{
		cutrank = maxrank;
	}

	double * sigma = (double *) malloc (maxrank * sizeof(double));
	double * tmp_factor = (double *) malloc (maxcarsize * sizeof(double));
	double * buffer = (double *) malloc(maxcarsize * sizeof(double));
	double * buffer2 = (double *) malloc(maxcarsize * sizeof(double));

	int * newranks = (int *) malloc ((this->dimensionality + 1) * sizeof(int));
	newranks[this->dimensionality] = this->ranks[this->dimensionality];
	newranks[0] = this->ranks[0];

	// <----- right_to_left - orthogonalization
	copy (this->ranks[dimensionality] * this->modes_sizes[dimensionality - 1] * this->ranks[dimensionality - 1], this->carriages[dimensionality - 1], 1, buffer2, 1);
    
	for (int i = dimensionality - 1; i > 0; i--)
	{
		newranks[i] = min(this->ranks[i], newranks[i + 1] * this->modes_sizes[i]);

		// now buffer2 is [newrank[i+1] x this->ranks[i] * this->modes_sizes[i] ]
		//      copy carriage_i to this->carriages[i]: [ newranks[i + 1] x ranks [i] x modes_sizes[i]]  -> [ranks[i] x modes_sizes[i] x newranks[i + 1]]
		transpose(buffer2, newranks[i + 1] * this->ranks[i], this->modes_sizes[i], this->carriages[i]);

		newranks[i] = min(newranks[i + 1] * this->modes_sizes[i], this->ranks[i]);

		// qr for this->carriages[i]: (modes_sizes[i] * newranks[i+1])  x  ranks[i]
		my_qr (this->modes_sizes[i] * newranks[i + 1], this->ranks[i], this->carriages[i], this->modes_sizes[i] * newranks[i + 1],  tmp_factor, newranks[i]);

		this->carriages[i] = (double *) realloc(this->carriages[i], this->modes_sizes[i] * newranks[i] * newranks[i + 1]* sizeof(double));
		//      now  this->carriages[i]  has sizes   modes_sizes[i] x newrank[i+1] x newrank[i]

		gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, newranks[i], this->ranks[i - 1] * this->modes_sizes[i-1], this->ranks[i], 1.0,
				tmp_factor,  newranks[i], this->carriages[i - 1], this->ranks[i], 0.0, buffer2, newranks[i]);
		// now buffer2 is [newrank[i] x this->ranks[i - 1] * this->modes_sizes[i-1] ]

		//       cout << "<--:  newrank [" << i << "] = " << newranks[i] << endl;
	}

	copy (this->dimensionality + 1, newranks, 1, this->ranks, 1);

    
	int zero = 0.0;
	copy(dimensionality + 1, &zero, 0, newranks, 1);
	newranks[0] = this->ranks[0];
	// carriage[0] (buffer2) - ranks[1] x ranks[0] x modes_sizes[0]
	// carriage[i] - modes_sizes[i] x ranks[i+1] x ranks[i]  for i > 0

	// -----> left_to_right - compression

//    cout << "myreducedorthofactors parameters: " << this->ranks[1] << " " << this->ranks[0] << " " << this->modes_sizes[0] << endl;
	myreducedorthofactors(eps, this->ranks[1], this->ranks[0] * this->modes_sizes[0], buffer2, tmp_factor, this->carriages[0], newranks[1], cutrank);
	// carriage[0] - newranks[1] x ranks[0] x  modes_sizes[0]
	// tmp_factor - ranks[1] x newranks[1]


	this->carriages[0] = (double *) realloc(this->carriages[0], newranks[1] * ranks[0] * modes_sizes[0] * sizeof(double));
	
    //    cout << "-->:  newrank [" << 1 << "] = " << newranks[1] << endl;
	for (int i = 1; i < this->dimensionality - 1; i++)
	{
		// carriage[i] - modes_sizes[i] x ranks[i+1] x ranks[i]
		// tmp_factor - ranks[i] x newranks[i]

		gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, modes_sizes[i] * this->ranks[i + 1], newranks[i], this->ranks[i], 1.0,
				this->carriages[i], modes_sizes[i] * this->ranks[i + 1], tmp_factor, this->ranks[i], 0.0, buffer, modes_sizes[i] * this->ranks[i+1]);
		// buffer - modes_sizes[i] x ranks[i+1] x newranks[i]

		transpose(buffer, this->modes_sizes[i], this->ranks[i + 1] * newranks[i], buffer2);
		// buffer2 - ranks[i+1] x newranks[i] x modes_sizes[i] 

		myreducedorthofactors(eps, this->ranks[i + 1],  newranks[i] * this->modes_sizes[i], buffer2, tmp_factor, this->carriages[i], newranks[i+1], cutrank);

		// carriages[i] -   newranks[i+1] x newranks[i] x modes_sizes[i]
		// tmp_factor  - ranks[i + 1] x newranks[i + 1]

		this->carriages[i] = (double *) realloc(this->carriages[i], newranks[i + 1] * newranks[i]  * modes_sizes[i] * sizeof(double));

		//cout << "-->:  newrank [" << i+1 << "] = " << newranks[i + 1] << endl;
	}
	newranks[dimensionality] = this->ranks[dimensionality];

	gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, modes_sizes[dimensionality-1] * this->ranks[dimensionality], newranks[dimensionality-1], this->ranks[dimensionality-1],
	1.0, this->carriages[dimensionality-1], modes_sizes[dimensionality-1] * this->ranks[dimensionality], tmp_factor, this->ranks[dimensionality-1], 0.0, buffer, modes_sizes[dimensionality - 1] * this->ranks[dimensionality]);

	// buffer - modes_sizes[d-1] x ranks[d] x newranks[d-1]
    
	transpose(buffer, this->modes_sizes[dimensionality - 1], this->ranks[dimensionality] * newranks[dimensionality - 1], this->carriages[dimensionality - 1]);

	this->carriages[dimensionality-1] = (double *) realloc(this->carriages[dimensionality - 1], newranks[dimensionality - 1] * modes_sizes[dimensionality - 1] * this->ranks[dimensionality] * sizeof(double));

	copy(this->dimensionality, newranks, 1, this->ranks, 1);


	free(buffer);
	free(buffer2);
	free(tmp_factor);
	free(sigma);
	free(newranks);

	return;
}


TTensorTrain & TTensorTrain::operator= (const TTensorTrain &tens)
{
	if (this == &tens) 
	{
		return *this;
	}
	free(ranks);
	free(modes_sizes);
	if (carriages != NULL) 
	{
		for (int i = 0; i < dimensionality; i++) 
		{
			free(carriages[i]);	
		}
		free(carriages);	
	}
	this->modes_sizes = NULL;
	this->ranks = NULL;
	this->carriages = NULL;

	this->dimensionality = tens.dimensionality;

	if (tens.modes_sizes == NULL)
	{
		return * this;
	}
	this->modes_sizes = (int *) malloc (dimensionality * sizeof(int));
	copy(this->dimensionality, tens.modes_sizes, 1, this->modes_sizes, 1);
   
	if (tens.ranks == NULL) 
	{
		return * this;
	}
	this->ranks = (int *) malloc ((dimensionality + 1) * sizeof(int));
	copy(this->dimensionality + 1, tens.ranks, 1, this->ranks, 1);
   
	this->carriages = (double **) malloc (this->dimensionality * sizeof (double *));
	for (int i = 0; i < this->dimensionality; ++i) 
	{
		this->carriages[i] = (double *) malloc(this->ranks[i] * this->modes_sizes[i] * this->ranks[i + 1] * sizeof(double));            
		copy(this->ranks[i] * this->modes_sizes[i] * this->ranks[i + 1], tens.carriages[i], 1, this->carriages[i], 1);
	}
	return *this;
}


TTensorTrain TTensorTrain :: operator+(const TTensorTrain &tensor)  const
{
	TTensorTrain sum = TTensorTrain();
	if (this->dimensionality != tensor.dimensionality) 
	{
		printf("Achtung! operator+ throws 1.\n");
		throw(1);
	}
	sum.dimensionality = tensor.dimensionality;

	sum.modes_sizes = (int *) malloc(sum.get_dimensionality() * sizeof(int));

	for (int i = 0; i < this->dimensionality; ++i) 
	{
		if (this->modes_sizes[i] != tensor.modes_sizes[i]) 
		{   
			printf("Achtung! Operator+ throws 2.\n");
			throw(2);
		}
		sum.modes_sizes[i] = this->modes_sizes[i];
	}

	sum.ranks = (int *) malloc((sum.get_dimensionality() + 1) * sizeof(int));
	sum.ranks [0] = 1;
	for (int i = 1; i < this->dimensionality; ++i) 
	{
		sum.ranks[i] = this->ranks[i] + tensor.ranks[i];
	}
	sum.ranks[this->dimensionality] = 1;

	sum.carriages = (double **) malloc(sum.dimensionality * sizeof(double*));
	double zero = 0.0;
	int size = sum.modes_sizes[0] * sum.ranks[1];
	sum.carriages[0] = (double *) malloc(size * sizeof(double));

	for (int j=0; j<sum.modes_sizes[0]; ++j)
	{
		copy(this->ranks[1], this->carriages[0] + j * this->ranks[1], 1, sum.carriages[0] + j * sum.ranks[1], 1);
		copy(tensor.ranks[1], tensor.carriages[0] + j * tensor.ranks[1], 1, sum.carriages[0] + j * sum.ranks[1] + this->ranks[1], 1);
	}

	for (int i = 1; i < dimensionality - 1; ++i)
	{
		size = sum.modes_sizes[i] * sum.ranks[i] * sum.ranks[i+1];
		sum.carriages[i] = (double *) malloc(size * sizeof(double));
		for (int j=0; j < sum.modes_sizes[i]; ++j)
		{
			for (int alpha = 0; alpha < this->ranks[i]; alpha++)
			{
				copy(this->ranks[i + 1], this->carriages[i] + j * this->ranks[i] * this->ranks[i + 1] + alpha * this->ranks[i + 1], 1, 
						sum.carriages[i] + j * sum.ranks[i] * sum.ranks[i + 1] + alpha * sum.ranks[i + 1], 1);

				copy(tensor.ranks[i + 1], &zero, 0, sum.carriages[i] + j * sum.ranks[i] * sum.ranks[i+1] + alpha * sum.ranks[i+1] + this->ranks[i+1], 1);

			}

			for (int alpha = this->ranks[i]; alpha < sum.ranks[i]; alpha++)
			{
				copy(this->ranks[i+1], &zero, 0, sum.carriages[i] + j * sum.ranks[i] * sum.ranks[i+1] + alpha * sum.ranks[i+1], 1);
				copy(tensor.ranks[i+1], tensor.carriages[i] + j * tensor.ranks[i] * tensor.ranks[i+1] + (alpha - this->ranks[i]) * tensor.ranks[i+1],
						1, sum.carriages[i] + j * sum.ranks[i] * sum.ranks[i+1] + alpha * sum.ranks[i+1] + this->ranks[i+1], 1);
			}
		}
	}

	size = sum.modes_sizes[sum.dimensionality - 1] * sum.ranks[sum.dimensionality - 1];
	sum.carriages[sum.dimensionality - 1] = (double *) malloc(size * sizeof(double));

	for (int j = 0; j < sum.modes_sizes[sum.dimensionality-1]; ++j)
	{
		copy(this->ranks[sum.dimensionality - 1], this->carriages[sum.dimensionality - 1] + j * this->ranks[sum.dimensionality - 1], 1, 
				sum.carriages[sum.dimensionality - 1] + j * sum.ranks[sum.dimensionality - 1], 1);            

		copy(tensor.ranks[sum.dimensionality - 1], tensor.carriages[sum.dimensionality - 1] + j * tensor.ranks[sum.dimensionality - 1], 1, 
				sum.carriages[sum.dimensionality - 1] + j * sum.ranks[sum.dimensionality - 1] + this->ranks[sum.dimensionality - 1], 1);
	}  
	//    std::cout << "95!" << std::endl;
	return sum;
}

TTensorTrain TTensorTrain :: operator-(const TTensorTrain &tensor)  const
{
    return * this + (tensor * (-1.0));
}

TTensorTrain TTensorTrain :: operator*(double alpha) const
{
    TTensorTrain new_tensor = TTensorTrain();
    new_tensor = (* this);
    int lastsize = new_tensor.modes_sizes[new_tensor.dimensionality - 1] * new_tensor.ranks[new_tensor.dimensionality-1] * new_tensor.ranks[new_tensor.dimensionality];
    scal (lastsize, alpha, new_tensor.carriages[new_tensor.dimensionality - 1], 1);
    return new_tensor;
}


TTensorTrain TTensorTrain :: operator *(const TTensorTrain &tensor) const
{
    int d = this->get_dimensionality();
    int *modes_sizes = (int *) malloc (d * sizeof(int));
    
    for (int i = 0; i < d; i++)
        modes_sizes[i] = this->get_mode_size(i);
    
    //TTensorTrain copy_of_this = &(this);
    MulTTensor multiplication_tensor(d, modes_sizes);
    multiplication_tensor.set_trains(this, &(tensor));
    
    TTensorTrainParameters parameters;
    parameters.maximal_iterations_number = 0;
    parameters.tolerance = 1e-6;
    
    TTensorTrain multiplication_tt;
    multiplication_tt.Approximate(&multiplication_tensor, parameters);
    return multiplication_tt;
}

TTensorTrain TTensorTrain :: elementwise_product(const TTensorTrain &tensor, const TTensorTrainParameters &parameters, const int if_compress, const int cutrank, const double eps_cut)
{
    int d = this->get_dimensionality();
    int *modes_sizes = (int *) malloc (d * sizeof(int));
    
    for (int i = 0; i < d; i++)
        modes_sizes[i] = this->get_mode_size(i);
    
    //TTensorTrain copy_of_this = &(this);
    MulTTensor multiplication_tensor(d, modes_sizes);
    multiplication_tensor.set_trains(this, &(tensor));
    
    TTensorTrain multiplication_tt;
    multiplication_tt.Approximate(&multiplication_tensor, parameters);
    
    
    if (if_compress)
        multiplication_tt.SVDCompress(eps_cut, cutrank);
    return multiplication_tt;
}
