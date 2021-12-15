#include "cross.h"

TVolume::TVolume(const int & s, const double & v, int * &i, int * &j)
{
    size = s;
    rows_positions = (int *) malloc(size*sizeof(int));
    columns_positions = (int *) malloc(size*sizeof(int));
    copy(size, i, 1, rows_positions, 1);
    copy(size, j, 1, columns_positions, 1);
    volume = v;
}


TVolume::TVolume(const int & s)
{
    size = s;
    rows_positions = (int *) calloc(size, sizeof(int));
    columns_positions = (int *) calloc(size, sizeof(int));
}

TVolume::~TVolume()
{
    free(rows_positions);
    free(columns_positions);
}

void TVolume::set(const double & v, int * i, int * j)
{
    copy(size, i, 1, rows_positions, 1);
    copy(size, j, 1, columns_positions, 1);
    volume = v;
}

double TVolume::get_volume() const
{
    return volume;
}

int TVolume::get_row_position() const
{
    return rows_positions[0];
}

int TVolume::get_column_position() const
{
    return columns_positions[0];
}

bool TVolume::operator>(const TVolume & vol) const
{
    return (abs(volume)>abs(vol.get_volume()))?true:false;
}

bool TVolume::operator==(const TVolume & vol) const
{
    return ((abs(volume) == abs(vol.get_volume())) && (rows_positions[0] == vol.get_row_position()) && (columns_positions[0] == vol.get_column_position()))?true:false;
}

bool TVolume::operator!=(const TVolume & vol) const
{
    return !((*this) == vol);
}

void TVolume::operator=(const TVolume & vol)
{
    volume = vol.get_volume();
    rows_positions[0] = vol.get_row_position();
    columns_positions[0] = vol.get_column_position();
}


