#include <algorithm>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <iomanip>
#include <cassert>
#include <utility>
#include "wrappers.h"
#include <ctime>
#include <cstdlib>
#include <limits>
#include <iterator>


using namespace std;
using namespace wrappers;

template<typename T>
std::vector<T> split(const std::string& line) {
    std::istringstream is(line);
    return std::vector<T>(std::istream_iterator<T>(is), std::istream_iterator<T>());
}

//solve first Smoluchowski integral
double L1(const int &N, const int &i, const double *n)
{
    double l1 = 0;
    for (int i1 = 0; i1 < i; i1++)
    {
        l1 += n[i1] * n[i - i1 - 1] * wrappers::K((i - i1 - 1), i1, 1);
    }
    return l1;
}

//solve second Smoluchowski integral
double L2(const int &N, const int &i, const double *n)
{
    double l2 = 0;
    for (int i1 = 0; i1 < N; i1++)
    {
        l2 += n[i1] * wrappers::K(i, i1, 1);
    }
    return l2;
}

//solve tridiagonal linear system
void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
    double m;
    for (int i = 1; i < n; i++)
    {
        m = a[i] / c[i-1];
        c[i] = c[i] - m * b[i-1];
        f[i] = f[i] - m * f[i-1];
    }
    x[n-1] = f[n-1]/c[n-1];

    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = ( f[i] - b[i] * x[i+1] ) / c[i];
    }
}

//use this to solve for case D = 0, but bugged, fix later
/*void only_advection_solver(double *&values, double const &dt, double *&coefs, int const &batch)
{
	double a, b;
	double ** res = new double * [size];
	for (int i = 0; i < size; i++) res[i] = new double [batch];
	for (int m = 0; m < batch; m++)
	{
		res[0][m] = values[0 + m * size];
		res[size - 1][m] = values[size - 1 + m * size];
	}
	for (int i = 1; i < size - 1; i++)
    {
		for (int m = 0; m < batch; m++)
        {
            b = values[i + 1 + m * size] - values[i - 1 + m * size];
            if (( values[i + 1 + m * size] - values[i + m * size] ) * ( values[i + m * size]-values[ i - 1 + m * size] ) > 0)
            {
                double first_val = min(abs(b) / ( 2.0 * dr ), 2.0 * abs(values[i + 1 + m * size] - values[i + m * size]) / dr);
                double second_val = 2.0 * abs(values[i + m * size] - values[i - 1 + m * size]) / dr;
                a = min(first_val, second_val) * ( b >= 0 ? 1.0 : -1.0 );
            }
            else
                a=0;
            res[i][m] = values[i + m * size] + dr / 2.0 * ( 1.0 - coefs[m] * dt / dr ) * a;
        }
    }
    for (int i = 1; i < size - 1; i++)
        for(int m = 0; m < batch; m++)
            values[i + m * size] = coefs[m] * ( ( res[i][m] - res[i - 1][m] ) / dr + sigma(dr * i, dr, PML, size) * values[i + m * size] );
    for (int i = 0; i < size; i++)
        delete [] res[i];
    delete [] res;

}*/

int main(int argc, char ** argv)
{
    ///*
    //res_full.txt

    double h = 1.0; 
    double dt = 0.01;
    double J = 1.0;
    int TIME = 20000;
    int MOD = 2000;
    int max_particle_size = 256;
    int max_x = 5000;
    double dx = 0.02;
    //*/

    /*
    //res_ballistic.txt

    double h = 1.0; 
    double dt = 0.001;
    double J = 1.0;
    int TIME = 50000;
    int MOD = 5000;
    int max_particle_size = 256;
    int max_x = 5000;
    double dx = 0.02;
    */

    /*
    //res_adv.txt
    
    double h = 1.0; 
    double dt = 0.1;
    double J = 1.0;
    int TIME = 5000; 
    int MOD = 500;
    int max_particle_size = 128;
    int max_x = 1000;
    double dx = 0.2;
    */

    //res_diff.txt
    /*
    double h = 1.0;
    double dt = 0.1;
    double J = 1.0;
    int TIME = 20000; 
    int MOD = 2000;
    int max_particle_size = 256;
    int max_x = 1000;
    double dx = 0.2;
    */

    double tolerance = 1e-3;    			
    TCross_Parallel_v1 crossed_kernel = default_crossed_kernel(tolerance, max_particle_size, h);		

    //FFTW START
    double R_value = crossed_kernel.get_rank();
    std::cout << "rank is " << R_value << std::endl;
    double V_value = crossed_kernel.get_columns_number();
    fftw_complex *ub = (fftw_complex *) fftw_malloc(R_value * V_value * sizeof(fftw_complex));
    fftw_complex *vb = (fftw_complex *) fftw_malloc(R_value * V_value * sizeof(fftw_complex));
    fftw_plan * plan_v = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    fftw_plan * plan_u = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    fftw_plan * plan_inverse = (fftw_plan *) fftw_malloc(R_value * sizeof(fftw_plan));
    for (int i = 0; i < R_value; i++)
    {
        plan_v[i] = fftw_plan_dft_1d(max_particle_size, vb+i*max_particle_size, vb+i*max_particle_size, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_u[i] = fftw_plan_dft_1d(max_particle_size, ub+i*max_particle_size, ub+i*max_particle_size, FFTW_FORWARD, FFTW_ESTIMATE);
        plan_inverse[i] = fftw_plan_dft_1d(max_particle_size, ub+i*max_particle_size, ub+i*max_particle_size, FFTW_BACKWARD, FFTW_ESTIMATE);	
    }
    //FFTW END

    double *n_k;
    double *L1_res_vec, *L2_res_vec;
    n_k = new double [max_particle_size];


    double *initial_layer = new double [max_particle_size * max_x];
    double *smoluch_operator = new double [max_particle_size * max_x];

    double *vel_coefs = new double [max_particle_size];
    double *dif_coefs = new double [max_particle_size];

    for (int m = 0; m < max_particle_size; m++)
    {
        //res_diff.txt
        /*
        vel_coefs[m] = 0.0;
        dif_coefs[m] = 1.0;
        */

        /*
        //res_adv.txt
        vel_coefs[m] = 1.0;
        dif_coefs[m] = 1.0;
        */

        ///*
        //res_full.txt or res_ballistic.txt
        vel_coefs[m] = 1.0*pow(m+1, 2./3.);
        dif_coefs[m] = 1.0*pow(m+1, -1./3.);
        //*/
        
    }
    double *coef_a = new double [max_x * max_particle_size];
    double *coef_b = new double [max_x * max_particle_size];
    double *coef_c = new double [max_x * max_particle_size];
    
    for (int m = 0; m < max_particle_size; m++)
    {
        for (int x = 0; x < max_x; x++)
        {
            int ind = x + m * max_x;
            coef_b[ind] = - dt * dif_coefs[m] * exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
            coef_c[ind] = 1.0 + dt * ( dif_coefs[m] * exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) + dif_coefs[m] * exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) ) / (dx*dx);
            coef_a[ind] = - dt * dif_coefs[m] * exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
        }

        coef_b[0 + m*max_x] = 1.0;

        coef_b[(max_x-1) + m*max_x] = 0.0;
        coef_a[0 + m*max_x] = 0.0;
        coef_a[(max_x-1) + m*max_x] = 0.0;

        coef_c[0 + m*max_x] = -1.0;

        coef_c[(max_x-1) + m*max_x] = 1.0;
    }

    clock_t start = clock();
    double duration;
    const auto stream_max = std::numeric_limits<std::streamsize>::max();
    ofstream output;
    ifstream input;
    if (argc !=3)
    {
        std::cout << "need output file name" << std::endl;
        return 2;
    }
    int append = atoi(argv[2]);
    if (append != 0){
        input.open(argv[1]);
        output.open(argv[1], std::ios_base::app);
        for (int i = 0; i < append*max_x; i++)
            input.ignore(stream_max, '\n');
        std::string line;
        int cc = 0;
        while (std::getline(input, line)){
            std::vector<double> vec = split<double>(line);
            for (int m = 0; m < max_particle_size; m++)
                initial_layer[m+max_particle_size*cc] = vec[m];
            cc++;
        }
        input.close();
    } else
        output.open(argv[1]);

    for (int t = 0; t < TIME; t++)
    {
        cout << t << endl;
        //this is calculating coagulation
        for (int x = 0; x < max_x; x++)
        {
            for (int m = 0; m < max_particle_size; m++)
            {
                int ind = m+x*max_particle_size;
                n_k[m] = initial_layer[ind];
                if (t%MOD==0) output << initial_layer[ind] << " ";
            }
            if (t%MOD==0) output << std::endl;
            
            
            //*****************************LOW RANK SOLUTION*************************************
            //clock_t st = clock();
            /*
            L2_res_vec = crossed_kernel.matvec(n_k);
            L1_res_vec = crossed_kernel.smol_conv_discrete(n_k, ub, vb, plan_v, plan_u, plan_inverse);

            #pragma omp parallel for
            for (int m = 0; m < max_particle_size; m++)
            {
                int ind = m+x*max_particle_size;
                if (m == 0 && x == 0)
                    smoluch_operator[ind] = (L1_res_vec[m]*0.5-n_k[m]*L2_res_vec[m])*dt+n_k[m];
                else
                    smoluch_operator[ind] = (L1_res_vec[m]*0.5-n_k[m]*L2_res_vec[m])*dt+n_k[m];
                if (smoluch_operator[ind] < 0.0) smoluch_operator[ind] = 0.0;
            }

            delete [] L2_res_vec;
            delete [] L1_res_vec;
            */
            //duration = (clock() - st) / (double)CLOCKS_PER_SEC;
            //cout << "low rank duration = " << duration << endl;
            //***********************************************************************************
            
            //*******************************DIRECT SOLUTION*************************************
            //st = clock();
			///*
            #pragma omp parallel for
            for (int m = 0; m < max_particle_size; m++)
            {
            	int ind = m+x*max_particle_size;
                smoluch_operator[ind] = ( L1(max_particle_size,m,n_k)*0.5-n_k[m]*L2(max_particle_size,m,n_k) )*dt+n_k[m];
                if (smoluch_operator[ind] < 0.0) smoluch_operator[ind] = 0.0;
            }
			//*/
            //duration = (clock() - st) / (double)CLOCKS_PER_SEC;
            //cout << "direct duration = " << duration << endl;
            //***********************************************************************************
        }

        //this is part calculating advection
        #pragma omp parallel for
        for (int m = 0; m < max_particle_size; m++)
        {
            double * rhs = new double [max_x];
            for (int x = 0; x < max_x; x++)
            {
                rhs[x] = smoluch_operator[m+x*max_particle_size];
                if (m==0 && x==0) rhs[x] = -J*dx*0.5;
                
                //remove this line if velocity is nonzero
                //else if (x==0) rhs[x] = 0.0;
            }
            double * output = new double [max_x];
            solveMatrix (max_x, &coef_a[m*max_x], &coef_c[m*max_x], &coef_b[m*max_x], rhs, output);
            
            for (int x = 0; x < max_x; x++)
            {
                initial_layer[m+x*max_particle_size] = output[x];
            }
            delete [] output;
            delete [] rhs;
        }

        //this is part of calculating advection
        #pragma omp parallel for
        for (int m = 0; m < max_particle_size; m++)
        {
            for (int x = 0; x < max_x; x++)
            {
                int ind = x + m * max_x;
                coef_b[ind] = - dt * dif_coefs[m] * exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
                coef_c[ind] = 1.0 + dt * ( dif_coefs[m] * exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) + dif_coefs[m] * exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) ) / (dx*dx);
                coef_a[ind] = - dt * dif_coefs[m] * exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
            }
            coef_b[0 + m*max_x] = 1.0;

            coef_b[(max_x-1) + m*max_x] = 0.0;
            coef_a[0 + m*max_x] = 0.0;
            coef_a[(max_x-1) + m*max_x] = 0.0;

            coef_c[0 + m*max_x] = -1.0;

            coef_c[(max_x-1) + m*max_x] = 1.0;
        }    
    }
    duration = (clock() - start) / (double)CLOCKS_PER_SEC;
    cout << "duration " << duration << endl;

    delete [] smoluch_operator;
    delete [] initial_layer;
    delete [] n_k;
    delete [] vel_coefs;
    delete [] dif_coefs;
    delete [] coef_a;
    delete [] coef_b;
    delete [] coef_c;

    for (int i = 0; i < R_value; i++)
    {
        fftw_destroy_plan(plan_v[i]);
        fftw_destroy_plan(plan_u[i]);
        fftw_destroy_plan(plan_inverse[i]);	
    }
    output.close();
    fftw_free(vb);
    fftw_free(ub);
    fftw_free(plan_u);
    fftw_free(plan_v);
    fftw_free(plan_inverse);
    return 0;
}

