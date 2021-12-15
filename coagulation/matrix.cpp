#include "matrix.h"
TMatrix::TMatrix(const int &m, const int &n)
{
    rows_number = m;
    columns_number = n;
}

TMatrix::~TMatrix()
{
    rows_number = 0;
    columns_number = 0;
}

int TMatrix::get_rows_number() const
{
    return rows_number;
}

int TMatrix::get_columns_number() const
{
    return columns_number;
}

TSubMatrix::TSubMatrix(const int & m, const int & n, int * r, int * c, TMatrix * matr):TMatrix(m,n)
{
    row_indeces = (int *) malloc(m*sizeof(int));
    copy(m, r, 1, row_indeces, 1);
    column_indeces = (int *) malloc(n*sizeof(int));
    copy(n, c, 1, column_indeces, 1);
    matrix = matr;
}

TSubMatrix::~TSubMatrix()
{
    free(row_indeces);
    free(column_indeces);
    matrix = NULL;
}

double TSubMatrix::value(const int & i, const int & j)
{
    return matrix->value(row_indeces[i],  column_indeces[j]);
}

TDifferenceMatrix::TDifferenceMatrix(TMatrix *m, TMatrix *s):TMatrix(m->get_rows_number(), m->get_columns_number())
{
    if ((m->get_rows_number() != s->get_rows_number()) || (m->get_columns_number() != s->get_columns_number()))
    {
        this->~TDifferenceMatrix();
        throw(trying_to_subtract_matrices_of_different_sizes);
    }
    else
    {
        Minuend_Matrix = m;
        Subtrahend_Matrix = s;
    }
}

double TDifferenceMatrix::value(const int &i, const int &j)
{
    return Minuend_Matrix->value(i,j) - Subtrahend_Matrix->value(i,j);
}

