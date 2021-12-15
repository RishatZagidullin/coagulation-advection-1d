#ifndef TENSOR_H
#define TENSOR_H
#include "matrix.h"
#include <cstdlib>
using namespace std;

class TTensor
{
    protected:
        int dimensionality;
        int *modes_sizes;
        TTensor();
    public:
        TTensor(const int &, int* &);
        ~TTensor();
        virtual double operator[](int *) = 0;
        int get_dimensionality() const;
        int get_mode_size(const int &) const;
        int get_Unfolding_Matrix_Rows_Number(const int &) const;
        int get_Unfolding_Matrix_Columns_Number(const int &) const;
};

class TUnfoldingSubMatrix:public TMatrix
{
    private:
        TTensor *tensor;
        int k;
        int **row_indeces, **column_indeces;
    public:
        double value(const int &, const int &);
        TUnfoldingSubMatrix(TTensor* , int , int , int , int **, int **);
        ~TUnfoldingSubMatrix();
};

#endif
