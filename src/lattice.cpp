#include <lattice.h>
#include <cmath>
#include <GSLpp/basic_math.h>
/*
Lattice::Lattice(const GSL::Vector& a, const GSL::Vector& b,
    const GSL::Vector& c)
 : lat(3,3), r_lat(3,3), scale(a.norm<double>()), volume(), bz_volume()
{
    lat[0] = a;
    lat[1] = b;
    lat[2] = c;
    lat *= 1./scale;

    volume = std::abs((a/scale).dot((b/scale).cross(c/scale)));
    r_lat[0] = 2*M_PI * (b/scale).cross(c/scale)/volume;
    r_lat[1] = 2*M_PI * (c/scale).cross(a/scale)/volume;
    r_lat[2] = 2*M_PI * (a/scale).cross(b/scale)/volume;

    bz_volume = std::abs(r_lat[0].dot(r_lat[1].cross(r_lat[2])));
}
*/
