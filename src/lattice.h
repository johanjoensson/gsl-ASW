#ifndef LATTICE_H
#define LATTICE_H

#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/vector.h"

/***************************************************************************//**
* A class used to represent lattices, e.g. Bravais lattices useful in describing
* crystals.\n
* Contains:\n
* __lat__ - GSL::Matrix constructed from the three lattice vectors (row major).
*\n
* __r_lat__ - GSL::Matrix constructed from the three reciprocal lattice vectors
* (row major).\n
* __scale__ - Scale factor to multiply each lattice vector with.\n
* __volume__ - Volume of the unit cell.\n
* __bz_volume__ - Volume of the first Brillouin zone.\n
*******************************************************************************/
class Lattice{

public:
    GSL::Matrix lat, r_lat;
    double scale, volume, bz_volume;

    Lattice();
    Lattice(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c);
    Lattice(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);
};
#endif // LATTICE_H
