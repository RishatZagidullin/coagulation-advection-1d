#ifdef FFTW
#include "convolution.h"
namespace convolution_helpers
{
	void ComplexScale_fftw(fftw_complex a, double s)
	{
		a[0] = s * a[0];
		a[1] = s * a[1];
	}

	void ComplexMul_fftw(fftw_complex a, const fftw_complex b)
	{
		fftw_complex c;
		c[0] = a[0] * b[0] - a[1] * b[1];
		c[1] = a[0] * b[1] + a[1] * b[0];
		a[0] = c[0];
		a[1] = c[1];
	}

	void ComplexPointwiseMulAndScale_fftw(fftw_complex * a, const fftw_complex * b, int size, double scale)
	{
		for (int i = 0; i < size; i ++)
		{
			ComplexMul_fftw(a[i], b[i]);
			ComplexScale_fftw(a[i], scale);
		}
	}
}
#endif
