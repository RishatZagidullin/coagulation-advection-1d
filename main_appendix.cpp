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
double L1(const int &N, const int &i, const double *n, double h)
{
    double l1 = 0;
    for (int i1 = 0; i1 < i; i1++)
    {
        double val = n[i1] * n[i-i1-1] * wrappers::K_appendix((i-i1-1), i1, h);
        //if ((i-i1-1) != i1)
        //    val *= pow(1.0 + gamma_q, 5./3.);
        l1 += val;
    }
    return l1;
}

//solve second Smoluchowski integral
double L2(const int &N, const int &i, const double *n, double h)
{
    double l2 = 0;
    for (int i1 = 0; i1 < N; i1++)
    {
        double val = n[i1] * wrappers::K_appendix(i, i1, h);
        //if (i != i1)
        //    val *= pow(1.0 + gamma_q, 5./3.);
        l2 += val;
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

void readData(double * data, int SIZE, std::string inFileName) {

    std::ifstream inFile;
    inFile.open(inFileName.c_str());

    if (inFile.is_open())
    {
        for (int i = 0; i < SIZE; i++)
        {
            inFile >> data[i];
            data[i] += 1.;
        }

        inFile.close();
    }
    else {
        cerr << "Can't find input file " << inFileName << endl;
    }
}

void fill_I(double * I, std::string name, int TIME)
{
    if (name == "smooth_23_july")
    {
        for (int i = 0; i < TIME/4; i++)
        {
            double value = ((double) i)/ ((double) (TIME/4-1));
            I[i] = 2.5*sin(3.141592654*value/2.+3.141592654/2.) + 0.54037054444 + 2.28;
        }
        for (int i = TIME/4; i < TIME/2; i++)
        {
            double value = ((double) i-TIME/4)/ ((double) (TIME/4-1));
            double value_x = value < 0.3 ? 0.3 : value;
            I[i] = (value_x+0.5)*sin(value*9.0 + 2.4) + 2.28;
        }
        for (int i = TIME/2; i < 3*TIME/4; i++)
        {
            double value = ((double) i-TIME/2)/ ((double) (TIME/4-1));
            I[i] = 2.5*sin(3.141592654*value/2.+3.141592654/2.) + 0.54037054444 + 2.28;
        }
        for (int i = 3*TIME/4; i < TIME; i++)
        {
            double value = ((double) i-3*TIME/4)/ ((double) (TIME/4-1));
            double value_x = value < 0.3 ? 0.3 : value;
            I[i] = (value_x+0.5)*sin(value*9.0 + 2.4) + 2.28;
        }
    }
    if (name == "direct_23_july")
    {
        readData(I, 800, "data/23july_resampled.txt");
        readData(I+800, TIME-800, "data/23july_resampled.txt");
    }
    if (name == "smooth_25_july")
    {
        for (int i = 0; i < TIME/4; i++)
        {
            double value = ((double) i)/ ((double) (TIME/2-1));
            I[i] = 1.*(sin(value*7*0.897597901+0.9)+2.);
        }
        for (int i = TIME/4; i < TIME/2; i++)
        {
            double value = ((double) i-TIME/4)/ ((double) (TIME/4-1));
            I[i] = 1.*(sin(value*7*0.897597901+0.9)+2.);
        }
        for (int i = TIME/2; i < 3*TIME/4; i++)
        {
            double value = ((double) i-TIME/2)/ ((double) (TIME/4-1));
            I[i] = 1.*(sin(value*7*0.897597901+0.9)+2.);
        }
        for (int i = 3*TIME/4; i < TIME; i++)
        {
            double value = ((double) i-3*TIME/4)/ ((double) (TIME/4-1));
            I[i] = 1.*(sin(value*7*0.897597901+0.9)+2.);
        }
    }
    if (name == "direct_25_july")
    {
        readData(I, TIME/2, "data/25july_resampled.txt");
        readData(I+TIME/2, TIME/2, "data/25july_resampled.txt");
    }
    if (name == "3_days")
    {
        readData(I, TIME, "data/3days_resampled.txt");
    }
}

int main(int argc, char ** argv)
{
    double dt = 0.02;
    int MOD = 1;
    int X_MOD=10;
    int max_x = 1;//200;
    double dx = 0.5;

    std::string filename{argv[4]};
    std::ifstream arguments;
    arguments.open(filename);
    std::vector<std::string> data;
    std::string line;
    int ii = 0;
    while(getline(arguments, line))
    {
        data.push_back(line);
        ii++;
    }
    
    arguments.close();
    std::string name = data[0];
    int TIME = std::stoi(data[1]);
    int max_particle_size = std::stoi(data[2]);
    double h = std::stod(data[3]);
    int num_smol = std::stoi(data[4]);
    //num_smol = 6;
    double dif_coef = std::stod(data[5]);
    std::cout << name << " " << TIME << " " << max_particle_size << " " << h << " " << num_smol << " " << dif_coef << "\n";
    double *I = new double[TIME];
    fill_I(I, name, TIME);

    std::ofstream check;
    check.open(argv[3]);

    double alpha = 1.;
    double beta = 0.0001;
    double C_oz = 0.0;
    int k_max = 40;


    double *n_k;
    double *L1_res_vec, *L2_res_vec;
    n_k = new double [max_particle_size];


    double *initial_layer = new double [max_particle_size * max_x];

    double *vel_coefs = new double [max_particle_size];
    double *dif_coefs = new double [max_particle_size];

    for (int m = 0; m < max_particle_size; m++)
    {
        vel_coefs[m] = 1.0*pow(h*(m+1), 2./3.);
        dif_coefs[m] = dif_coef*pow(h*(m+1), -1./3.);
        n_k[m] = 0.0;
        for (int x = 0; x < max_x; x++)
        {
            initial_layer[m+x*max_particle_size] = 0.0;
        }
        
    }
    double *coef_a = new double [max_x * max_particle_size];
    double *coef_b = new double [max_x * max_particle_size];
    double *coef_c = new double [max_x * max_particle_size];
    
    for (int m = 0; m < max_particle_size; m++)
    {
        for (int x = 0; x < max_x; x++)
        {
            int ind = x + m * max_x;
            coef_b[ind] = - dt*dif_coefs[m]*exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
            coef_c[ind] = 1.0 + dt*( dif_coefs[m]*exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) +
                                     dif_coefs[m]*exp(-vel_coefs[m]/dif_coefs[m]*dx/2.0) ) / (dx*dx);
            coef_a[ind] = - dt*dif_coefs[m]*exp(vel_coefs[m]/dif_coefs[m]*dx/2.0) / (dx*dx);
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
    int append = atoi(argv[2]);
    if (append != 0){
        input.open(argv[1]);
        output.open(argv[1], std::ios_base::app);
        for (int i = 0; i < append*max_x; i++)
            input.ignore(stream_max, '\n');
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

    int count = -1;
    for (int t = 0; t < TIME; t++)
    {
        //if (t%MOD==0) cout << t << " ";
        for (int numm = 0; numm < num_smol; numm ++)
        {
            for (int x = 0; x < max_x; x++)
            {
                for (int m = 0; m < max_particle_size; m++)
                {
                    int ind = m+x*max_particle_size;
                    n_k[m] = initial_layer[ind];
                }
            
                //*****************DIRECT SOLUTION*****************
                #pragma omp parallel for
                for (int m = 0; m < max_particle_size; m++)
                {
                    int ind = m+x*max_particle_size;
                    initial_layer[ind] = ( L1(max_particle_size,m,n_k, h)*0.5 - 
                                           n_k[m]*L2(max_particle_size,m,n_k, h) )*dt+n_k[m];
                    if (initial_layer[ind] < 0.0) initial_layer[ind] = 0.0;
                }
                for (int m = 0; m < max_particle_size; m++)
                {
                    int ind = m+x*max_particle_size;
                    if (t%MOD==0 && x%X_MOD==0 && numm==num_smol-1) output << initial_layer[ind] << " ";
                }
                if (t%MOD==0 && x%X_MOD==0 && numm ==num_smol-1) output << std::endl;
                //*************************************************
            }
        }

        //modelling Baikal data
        double member = beta*C_oz;
        C_oz = C_oz + dt*(I[t]-alpha*C_oz-member);
        if (C_oz < 0) C_oz = 0.0;
        
        for (int x = 0; x < max_x; x++)
        {
            if (t%MOD==0 && x==0) check << C_oz << " ";
            double sums = 0.0;
            for (int m = 0; m < k_max; m++)
            {
                sums +=wrappers::K_appendix(0, m, h) *
                            initial_layer[m+x*max_particle_size] * 
                            initial_layer[0+x*max_particle_size];
            }
            initial_layer[0+x*max_particle_size] = 
                              initial_layer[0+x*max_particle_size] + 
                                          dt*(member - sums);
            //std::cout << initial_layer[0+x*max_particle_size] << " " << sums << " " << member << "\n";
            if (initial_layer[0+x*max_particle_size] < 0.0) initial_layer[0+x*max_particle_size] = 0.0;
        }
        if (t%MOD==0) check << "\n";
        //************************
    }

    delete [] initial_layer;
    delete [] n_k;
    delete [] vel_coefs;
    delete [] dif_coefs;
    delete [] coef_a;
    delete [] coef_b;
    delete [] coef_c;
    delete [] I;

    output.close();
    check.close();
    return 0;
}
