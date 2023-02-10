#include <algorithm>
#include <memory>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <utility>
#include "wrappers.h"
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
double L1(const int &N, const int &i, const double *n, double q)
{
    double l1 = 0;
    for (int i1 = 0; i1 < i; i1++)
    {
        l1 += n[i1] * n[i - i1 - 1] * wrappers::K((i - i1 - 1), i1, 1)*(1.0 + q);
    }
    return l1;
}

//solve second Smoluchowski integral
double L2(const int &N, const int &i, const double *n, double q)
{
    double l2 = 0;
    for (int i1 = 0; i1 < N; i1++)
    {
        l2 += n[i1] * wrappers::K(i, i1, 1)*(1.0+q);
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

void readData(double * data, int SIZE) {


    std::string inFileName = "data/2_resampled.txt";
    std::ifstream inFile;
    inFile.open(inFileName.c_str());

    if (inFile.is_open())
    {
        for (int i = 0; i < SIZE; i++)
        {
            inFile >> data[i];
            //std::cout << data[i] << " ";
        }

        inFile.close(); // CLose input file
    }
    else { //Error message
        cerr << "Can't find input file " << inFileName << endl;
    }
}

int main(int argc, char ** argv)
{
    ///*
    //res_full.txt

    double h = 1.0; 
    double dt = 0.001;//0.002;//0.005
    double J = 1.0;
    int TIME = 20000;//10000;//5000?
    int MOD = 20;
    int max_particle_size = 128;//64;//32
    int max_x = 500;//1000;//500
    double dx = 0.1;//0.05;//0.1
    //*/

    /*
    //res_ballistic.txt

    double h = 1.0; 
    double dt = 0.001;
    double J = 1.0;
    int TIME = 50001;
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
    int TIME = 5001; 
    int MOD = 500;
    int max_particle_size = 128;
    int max_x = 1000;
    double dx = 0.2;
    */

    //res_diff.txt
    /*
    double h = 1.0;
    double dt = 0.2;
    double J = 1.0;
    int TIME = 2001; 
    int MOD = 200;
    int max_particle_size = 256;
    int max_x = 500;
    double dx = 0.4;
    */

    //params for modelling Baikal data
    double *I = new double[TIME];
    //readData(I, TIME);
    std::ofstream check;
    check.open(argv[3]);
    for (int i = 0; i < TIME; i++)
    {
        //I[i] = I[i]/11.9;
        I[i] = (sin(i*0.00035+0.5)+1.);
        //I[i] = 0.5*exp(-0.000005*(i-TIME/2-TIME/17)*(i-TIME/2-TIME/17));
        //I[i] += 0.333*exp(-0.0000005*(i-TIME/2+TIME/25)*(i-TIME/2+TIME/25));
    }
    //for (int i = TIME/2; i < TIME; i++)
    //{
    //    I[i] = 0.5*exp(-0.00001*(i-TIME/2-TIME/6)*(i-TIME/2-TIME/6));
    //}
    //64
    double alpha = 5.0;
    double beta = 3.0;
    double gamma = 10.0;
    double q_max = 15.0;
    double ksi = 35.0;
    
    //128
    //double alpha = 1.0;
    //double beta = 0.5;
    //double gamma = 10.0;
    //double q_max = 15.0;
    //double ksi = 35.0;


    double *ozone = new double[max_x];
    //double *ozone = new double[TIME];
    //readData(ozone, TIME);
    
    double *q = new double[max_x];
    //********************************

    double tolerance = 1e-3;    			
    TCross_Parallel_v1 crossed_kernel = default_crossed_kernel(tolerance, max_particle_size, h);		

    //FFTW START
    double R_value = crossed_kernel.get_rank();
    std::cout << "rank is " << R_value << std::endl;
    std::cout << alpha << " " << beta << " " << gamma << "\n";
    std::cout << argv[3] << std::endl;
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
        /*
        //res_diff.txt
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
        vel_coefs[m] = 1.*pow(m+1, 2./3.);
        dif_coefs[m] = 1.*pow(m+1, -1./3.);
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

    const auto stream_max = std::numeric_limits<std::streamsize>::max();
    ofstream output;
    ifstream input;
    //if (argc !=3)
    //{
    //    std::cout << "need output file name" << std::endl;
    //    return 2;
    //}
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
        if (t%MOD==0) cout << t << "/" << TIME-1 << endl;
        //this is calculating coagulation
        for (int x = 0; x < max_x; x++)
        {
            for (int m = 0; m < max_particle_size; m++)
            {
                int ind = m+x*max_particle_size;
                n_k[m] = initial_layer[ind];
                if (t%MOD==0 && x%25==0) output << initial_layer[ind] << " ";
            }
            if (t%MOD==0 && x%25==0) output << std::endl;
            
            
            //**************LOW RANK SOLUTION******************
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
            //*************************************************
            
            //*****************DIRECT SOLUTION*****************
            ///*
            #pragma omp parallel for
            for (int m = 0; m < max_particle_size; m++)
            {
            	int ind = m+x*max_particle_size;
                smoluch_operator[ind] = ( L1(max_particle_size,m,n_k,q[x])*0.5-n_k[m]*L2(max_particle_size,m,n_k,q[x]) )*dt+n_k[m];
                if (smoluch_operator[ind] < 0.0) smoluch_operator[ind] = 0.0;
            }
            //*/
            //*************************************************
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
                
                //uncomment this line if res_diff.txt else comment it
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

        //modelling Baikal data
        for(int gg = 0; gg < 1; gg++)
        {
        for (int x = 0; x < max_x; x++)
        {
            double sums = 0.0;
            for (int m = 0; m < max_particle_size; m++)
            {
                sums += initial_layer[m+x*max_particle_size]*pow(m, 2./3.);
            }
            ozone[x] = ozone[x]+dt*(I[t]-alpha*ozone[x]-beta*sums*ozone[x]);
            if (t%MOD==0) check << ozone[x] << " ";
            if (ozone[x] < 0.0)
                ozone[x] = 0.0;
            q[x] = q[x]+dt*(gamma*ozone[x]*(q_max-q[x])*(q[x]>=q_max?0.:q_max)-ksi*q[x]);
            //q[x] = q[x]+dt*(gamma*ozone[t]/11.9*(q_max-q[x])*(q[x]>=q_max?0.:q_max)-ksi*q[x]);
            if (q[x] < 0.0)
                q[x] = 0.0;
        }
        }
        if (t%MOD==0) check << "\n";
        //************************

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

    delete [] smoluch_operator;
    delete [] initial_layer;
    delete [] n_k;
    delete [] vel_coefs;
    delete [] dif_coefs;
    delete [] coef_a;
    delete [] coef_b;
    delete [] coef_c;
    delete [] q;
    delete [] ozone;
    delete [] I;

    for (int i = 0; i < R_value; i++)
    {
        fftw_destroy_plan(plan_v[i]);
        fftw_destroy_plan(plan_u[i]);
        fftw_destroy_plan(plan_inverse[i]);	
    }
    output.close();
    check.close();
    fftw_free(vb);
    fftw_free(ub);
    fftw_free(plan_u);
    fftw_free(plan_v);
    fftw_free(plan_inverse);
    return 0;
}

