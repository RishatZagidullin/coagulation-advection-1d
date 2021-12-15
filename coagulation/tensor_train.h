#ifndef TENSOR_TRAIN_H
#define TENSOR_TRAIN_H
#include "tensor.h"
#include <cstdlib>
#include <cstdio>
#include "parallel_cross_omp.h"
#include <cmath>
#include <algorithm>
#include "blas.h"

struct TTensorTrainParameters
{
    public:
        double tolerance;
        int maximal_iterations_number, stop_rank, cross_max_iterations;
};

class TTensorTrain : public TTensor
{
    protected:
        int *ranks;
        double **carriages;
    public:
        double operator[](int *);
 
	double operator[](const int &);

	TTensorTrain operator*(const TTensorTrain &) const;
        TTensorTrain operator*(double) const;
        friend TTensorTrain operator*(double alpha, const TTensorTrain &tt)
        {
            return tt * alpha;
        }
	TTensorTrain operator+(const TTensorTrain &) const;
	TTensorTrain operator-(const TTensorTrain &) const;

	TTensorTrain & operator=(const TTensorTrain &);
	void Approximate(TTensor *tensor, const TTensorTrainParameters &);

        TTensorTrain elementwise_product(const TTensorTrain &, const TTensorTrainParameters &, const int if_compress,
                                         const int cutrank, const double eps_cut);

        int get_rank(int &);
        void save_to_file(char *);
        void load_from_file(char *);
        //void Orthogonalize();
        //void Compress(const double &);
        double nrm2();
        double dot(const TTensorTrain &);
        TTensorTrain();
        TTensorTrain(const TTensorTrain &, const  double &);
        TTensorTrain(const TTensorTrain &tens, const  double *c);
        ~TTensorTrain();
        void SVDCompress (const double = 0.0, const int = 0);
};

class MulTTensor: public TTensor{

    private:
        TTensorTrain *tt_one, *tt_two;
    public:
        double operator [] (int *);
        MulTTensor();
        MulTTensor(const int &d, int *&m);
        void set_trains( const TTensorTrain *,  const TTensorTrain *);

};

#endif

