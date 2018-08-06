#ifndef LATTICE_H
#define LATTICE_H

#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/vector.h"

class Lattice{

public:
    GSL::Matrix lat, r_lat;
    double scale, volume, bz_volume;

    Lattice();
    Lattice(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c);
    Lattice(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);
};
#endif // LATTICE_H
