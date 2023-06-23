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
        double u1=pow( (u + 1.0)*h , 1.0/3.0);
        double v1=pow( (v + 1.0)*h , 1.0/3.0);
        double result = 0.0;
        if (u==v)
            result = (u1+v1)*(1./u1+1./v1);
        else
            result = pow((u1+v1), 2./3.)*((u+1)*h+(v+1)*h) / fabs(u1*u1-v1*v1) / pow((u+1)*h*(v+1)*h, 5./9.);
            //result = pow((u1+v1), 2)*fabs(u1*u1-v1*v1);
        
        //kernel for res_diff.txt and res_adv.txt
        //double result = 2.0;
        return 1.*result;
    }
}

