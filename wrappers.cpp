#include "wrappers.h"
#include <cmath>


namespace wrappers
{
    double K(const int & u, const int &v, const double h)
    {
        //simplified kernel for res_full.txt
        //double u1=pow( (u + 1.0) , 2.0/3.0);
        //double v1=pow( (v + 1.0) , 2.0/3.0);
        //double result = u1*v1;

        //ballistic kernel for res_ballistic.txt
        //double u1=pow( (u + 1.0) , 1.0/3.0);
        //double v1=pow( (v + 1.0) , 1.0/3.0);
        //double result = 0.0;
        //if (u==v)
        //    result = (u1+v1)*(1./u1+1./v1);
        //else
        //    result = pow((u1+v1), 2)*fabs(u1*u1-v1*v1);
        
        //kernel for res_diff.txt and res_adv.txt
        double result = 2.0;
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

