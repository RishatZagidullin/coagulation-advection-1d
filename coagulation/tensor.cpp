#include "tensor.h"


TTensor::TTensor(const int &d, int* &m)
{
    dimensionality = d;
    modes_sizes = m;
    m = NULL;
}

TTensor::TTensor()
{
    dimensionality = 0 ;
    modes_sizes = NULL;
}

TTensor::~TTensor()
{
    dimensionality = 0 ;
    free(modes_sizes);
    modes_sizes = NULL;
}

int TTensor::get_dimensionality() const
{
    return dimensionality;
}

int TTensor::get_mode_size(const int &i) const
{
    if (i < dimensionality)
    {
        try
        {
            return modes_sizes[i];
        }
        catch(...)
        {
            return 0 ;
        }
    }
    else return 0 ;
}

int TTensor::get_Unfolding_Matrix_Rows_Number(const int &i) const
{
    if (i > dimensionality) return 0 ;
    int rows_number = 1 , j;
    for (j = 0 ; j < i; j++) rows_number *= modes_sizes[j];
    return rows_number;
}

int TTensor::get_Unfolding_Matrix_Columns_Number(const int &i) const
{
    if (i > dimensionality) return 0 ;
    int columns_number = 1 , j;
    for (j = i; j < dimensionality; j++) columns_number *= modes_sizes[j];
    return columns_number;
}

TUnfoldingSubMatrix::TUnfoldingSubMatrix(TTensor* tens, int k1, int m, int n, int ** r, int ** c):TMatrix(m, n) 
{
    this->tensor = tens;
    this->k = k1;
    row_indeces = r;
    column_indeces = c;
}

TUnfoldingSubMatrix::~TUnfoldingSubMatrix()
{
    tensor = NULL;
    row_indeces = NULL;
    column_indeces = NULL;
}


double TUnfoldingSubMatrix::value(const int & i, const int & j) 
{
    int index[this->tensor->get_dimensionality()];
    for (int o = 0 ; o < k; o++) index[o] = row_indeces[i][o];
    for (int o = k; o < this->tensor->get_dimensionality(); o++) index[o] = column_indeces[j][this->tensor->get_dimensionality() - o - 1 ];
    return (*this->tensor)[index];
}

