#include "wrappers.h"
#include <cmath>


namespace wrappers
{
    double K(const int & u, const int &v, const double h)
    {
        //ballistic kernel for res_ballistic.txt
        double u1=pow( (u + 1.0)*h , 1.0/3.0);
        double v1=pow( (v + 1.0)*h , 1.0/3.0);
        double result = 0.0;
        result = (u1+v1)*(u1+v1)*(1./(u+1)*h+1./(v+1)*h);

        //simple kernel
        //double result = 2.;
        return result;
    }
    //******************************************************************************************
    //******************************************************************************************
    TCross_Parallel_v1 default_crossed_kernel(const double & tolerance, const int & size, const double & dm)
    {
        TCross_Parallel_v1_Parameters parameters;
        parameters.tolerance = tolerance;
        parameters.maximal_iterations_number = 0;
        TKernel kernel(size, size);
        kernel.h = dm;
        TCross_Parallel_v1 crossed_kernel;
        crossed_kernel.Approximate(&kernel, parameters);
        return crossed_kernel;
    }
}

