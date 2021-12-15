#pragma once
#include <stdio.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <fstream>
#include "coagulation/tensor_train.h"
#include <vector>

using namespace std;
namespace wrappers
{

	double K(const int & u, const int &v, const double h);

  class TKernel: public TMatrix
	{
		public:
               	double h;
               	TKernel (const int &m, const int &n) : TMatrix(m , n) {}
		double value (const int &i, const int &j) override {return h*K(i,j,h);}
	};

	TCross_Parallel_v1 default_crossed_kernel(const double & tolerance, const int & size, const double & dm);
}

