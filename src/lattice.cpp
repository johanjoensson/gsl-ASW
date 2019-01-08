#include "lattice.h"
#include "GSLpp/basic_math.h"

Lattice::Lattice()
 : lat(), r_lat(), scale(), volume(), bz_volume()
{}

Lattice::Lattice(const GSL::Vector& a, const GSL::Vector& b,
    const GSL::Vector& c)
 : lat(3,3), r_lat(3,3), scale(a.norm()), volume(0), bz_volume()
{
    lat[0] = a;
    lat[1] = b;
    lat[2] = c;
    lat *= 1./scale;

    volume = std::abs(GSL::dot(lat[0], GSL::cross(lat[1], lat[2])));
    r_lat[0] = 2*M_PI * GSL::cross(b/scale, c/scale)/volume;
    r_lat[1] = 2*M_PI * GSL::cross(c/scale, a/scale)/volume;
    r_lat[2] = 2*M_PI * GSL::cross(a/scale, b/scale)/volume;

    bz_volume = std::abs(GSL::dot(r_lat[0], GSL::cross(r_lat[1], r_lat[2])));
}
