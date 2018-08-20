#include "lattice.h"
#include "../../GSL-lib/src/basic_math.h"

Lattice::Lattice()
 : lat(), r_lat(), scale(), volume(), bz_volume()
{}

Lattice::Lattice(GSL::Vector& a, GSL::Vector& b, GSL::Vector& c)
 : lat(3,3), r_lat(3,3), scale(a.norm())
{
    lat[0] = a/scale;
    lat[1] = b/scale;
    lat[2] = c/scale;

    volume = abs(GSL::dot(lat[0], GSL::cross(lat[1], lat[2])));
    r_lat[0] = 2*M_PI * GSL::cross(b/scale, c/scale)/volume;
    r_lat[1] = 2*M_PI * GSL::cross(c/scale, a/scale)/volume;
    r_lat[2] = 2*M_PI * GSL::cross(a/scale, b/scale)/volume;

    bz_volume = abs(GSL::dot(r_lat[0], GSL::cross(r_lat[1], r_lat[2])));
}

Lattice::Lattice(const GSL::Vector& a, const GSL::Vector& b,
    const GSL::Vector& c)
 : lat(3,3), r_lat(3,3), scale(a.norm())
{
    lat[0] = a/scale;
    lat[1] = b/scale;
    lat[2] = c/scale;

    volume = abs(GSL::dot(lat[0], GSL::cross(lat[1], lat[2])));
    r_lat[0] = 2*M_PI * GSL::cross(b/scale, c/scale)/volume;
    r_lat[1] = 2*M_PI * GSL::cross(c/scale, a/scale)/volume;
    r_lat[2] = 2*M_PI * GSL::cross(a/scale, b/scale)/volume;

    bz_volume = abs(GSL::dot(r_lat[0], GSL::cross(r_lat[1], r_lat[2])));
}
