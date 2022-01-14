#ifndef LATTICE_H
#define LATTICE_H

#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"
#include "GSLpp/linalg.h"

/*
class Lattice{

public:
    GSL::Matrix lat, r_lat;
    double scale, volume, bz_volume;

    Lattice(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);
};
*/

template<size_t dim>
class Lattice_t
{
private:
    GSL::Matrix lat_m, recip_lat_m;
    double scale_m;
public:

    Lattice_t(GSL::Matrix&& m)
     : lat_m(std::move(m)), recip_lat_m(lat_m.num_rows(), lat_m.num_columns()), scale_m(GSL::norm(lat_m[0]))
    {
        this->lat_m *= 1./scale_m;
        this->recip_lat_m.copy( 2*M_PI*GSL::lu_inverse(lat_m) );
    }

    double scale() const {return scale_m;}
    GSL::Matrix lat() const {return scale_m*lat_m.cview();}
    GSL::Matrix recip_lat() const {return 1/scale_m*recip_lat_m.cview();}
};
#endif // LATTICE_H
