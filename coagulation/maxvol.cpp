#include "maxvol.h"
int maxvol(int n, int r, const double *factor, int *start_index, double eps, double *coef)
{
	int size = n*r;
	int i, j;
	int *tmp_ipiv;
	double *tmp_matrix, *tmp_row, *tmp_column, *pointer1, *pointer2;
	double alpha = 1.0, maxf;
	tmp_matrix = (double *)malloc(r*r*sizeof(double));
	tmp_ipiv = (int *)malloc(r*sizeof(int));
	for (i = 0; i < r; i++)
	{
		copy(n, factor+i*n, 1, coef+i, r);
	}
	for (i = 0; i < r; i++)
	{
		pointer2 = tmp_matrix+i*r;
		copy(r, factor+start_index[i], n, pointer2, 1);
	}
	gesv(LAPACK_COL_MAJOR, r, n, tmp_matrix, r, tmp_ipiv, coef, r);
	free(tmp_ipiv);
	free(tmp_matrix);
	tmp_row = (double *) malloc(n*sizeof(double));
	tmp_column = (double *) malloc(r*sizeof(double));
	maxf = eps+2;
	while (maxf > eps + 1.0l)
	{
		i = iamax(size, coef, 1);
		maxf = coef[i];
		j = i/r;
		i -= j*r;
		if (fabs(maxf) > eps + 1.0l)
		{
			copy(n, coef+i, r, tmp_row, 1);
			copy(r, coef+j*r, 1, tmp_column, 1);
			tmp_column[i] -= 1.0;
			start_index[i] = j;
			alpha = -1.0/maxf;
			ger(CblasColMajor, r, n, alpha, tmp_column, 1, tmp_row, 1, coef, r);
		}
	}
	free(tmp_row);
	free(tmp_column);
	return 0;
}
