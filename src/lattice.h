#ifndef LATTICE_H
#define LATTICE_H

#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"
#include "GSLpp/linalg.h"

class Lattice{

public:
    GSL::Matrix lat, r_lat;
    double scale, volume, bz_volume;

    // Lattice() : lat{}, r_lat{}, scale{}, volume{}, bz_volume{} {};
    Lattice(const GSL::Vector& a, const GSL::Vector& b, const GSL::Vector& c);
};


template<size_t dim>
class Lattice_t
{
private:
    GSL::Matrix lat_m, recip_lat_m;
    double scale_m;
public:
    // Lattice_t() : lat_m{}, recip_lat_m{}, scale_m{} {};

    Lattice_t(const std::initializer_list<GSL::Vector>& l)
     : lat_m(l), recip_lat_m(), scale_m(l.begin()->norm<double>())
    {
        lat_m *= 1./scale_m;
        recip_lat_m = 2*M_PI*GSL::lu_inverse(lat_m);
    }

    double scale() const {return scale_m;}
    GSL::Matrix lat() const {return scale_m*lat_m;}
    GSL::Matrix recip_lat() const {return 1/scale_m*recip_lat_m;}
};
#endif // LATTICE_H
