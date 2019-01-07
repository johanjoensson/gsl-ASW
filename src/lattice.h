#ifndef LATTICE_H
#define LATTICE_H

#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"

class Lattice{

public:
    GSL::Matrix lat, r_lat;
    double scale, volume, bz_volume;

    Lattice();
    Lattice(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);
};
#endif // LATTICE_H
