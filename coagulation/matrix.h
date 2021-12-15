#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
#include "errors.h"
#include "blas.h"

class TMatrix
{
    protected:
        int rows_number;
        int columns_number;
        TMatrix(const int &, const int &);
        ~TMatrix();
    public:
        virtual double value(const int &, const int &) = 0;
        int get_rows_number() const;
        int get_columns_number() const;
};

class TSubMatrix: public TMatrix
{
    private:
        int *row_indeces, *column_indeces;
        TMatrix *matrix;
    public:
        double value(const int &, const int &);
        TSubMatrix(const int &, const int &, int* , int* , TMatrix* matr);
        ~TSubMatrix();
};

class TDifferenceMatrix: public TMatrix
{
    private:
        TMatrix *Minuend_Matrix, *Subtrahend_Matrix;
    public:
        TDifferenceMatrix(TMatrix *, TMatrix *);
        double value(const int &, const int &);
};

#endif
