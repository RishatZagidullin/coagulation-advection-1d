#ifndef CROSS_H
#define CROSS_H
#include "matrix.h"
#include "blas.h"
#include <cmath>
using namespace std;

class TVolume{
    private:
        double volume;
        int *rows_positions, *columns_positions;
        int size;
    public:
        TVolume(const int &, const double &, int *&, int *&);
        TVolume(const TVolume &);
        TVolume(const int &);
        void set(const double &, int *, int *);
        double get_volume() const;
        int get_row_position() const;
        int get_column_position() const;
        bool operator>(const TVolume &) const;
        bool operator==(const TVolume &) const;
        void operator=(const TVolume &);
        bool operator!=(const TVolume &) const;
        ~TVolume();
};

template <class TCross_Work_Data, class TCross_Parameters>
class TCross_Base:public TMatrix
{
    protected:
        int rank;
        virtual void Prepare_Data(TCross_Work_Data &, const TCross_Parameters &) = 0;
        virtual void Search_Max_Volume(TCross_Work_Data &) = 0;
        virtual bool Stopping_Criteria(TCross_Work_Data &) = 0;
        virtual void Update_Cross(TCross_Work_Data &) = 0;
    public:
        void Approximate(TMatrix * matrix, const TCross_Parameters &parameters)
        {
            this->rows_number = matrix->get_rows_number();
            this->columns_number = matrix->get_columns_number();
            TCross_Work_Data WorkData(matrix, this);
            Prepare_Data(WorkData, parameters);
            while (true)
            {
                Search_Max_Volume(WorkData);
                bool b = Stopping_Criteria(WorkData);
                if (!b)
                {
                    Update_Cross(WorkData);
                }
                else
                {
                    break;
                }
            }
        }

        TCross_Base():TMatrix(0,0)
        {
        }

        int get_rank()
        {
            return rank;
        }
};

class TMax_Volumes_Array{
    private:
        TVolume *volumes;
        int number_of_volumes;
    public:
        int get_number_of_volumes();
        const TVolume & get_volume(const int &);
        const TVolume & operator[](const int &);
        void insert(const TVolume &);
};
#endif
