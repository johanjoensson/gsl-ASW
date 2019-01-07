#include "k-mesh.h"
#include <iostream>
#include <cmath>

K_mesh::K_mesh()
 : r_lattice(), k_points()
{}

K_mesh::K_mesh(const GSL::Matrix& r_lat)
 : r_lattice(r_lat), k_points()
{}

void K_mesh::generate_mesh(const size_t nx, const size_t ny, const size_t nz)
{
    const GSL::Vector Gx = r_lattice[0], Gy = r_lattice[1], Gz = r_lattice[2];
    GSL::Vector tmp(3);
    double nxd = static_cast<double>(nx);
    double nyd = static_cast<double>(ny);
    double nzd = static_cast<double>(nz);
    double twonxd = static_cast<double>(2*nx);
    double twonyd = static_cast<double>(2*ny);
    double twonzd = static_cast<double>(2*nz);

    for(size_t n1 = 0; n1 < nx; n1++){
        for(size_t n2 = 0; n2 < ny; n2++){
            for(size_t n3 = 0; n3 < nz; n3++){
                tmp = (static_cast<double>(2*n1) - nxd - 1)/twonxd*Gx;
                tmp += (static_cast<double>(2*n2) - nyd - 1)/twonyd*Gy;
                tmp += (static_cast<double>(2*n3) - nzd - 1)/twonzd*Gz;
                k_points.push_back({tmp[0], tmp[1], tmp[2]});
            }
        }
    }
}

void K_mesh::generate_mesh(const double r_max)
{
    const GSL::Vector Gx = r_lattice[0], Gy = r_lattice[1], Gz = r_lattice[2];
    size_t nx = static_cast<size_t>(std::max(r_max/Gx.norm(), 1.));
    size_t ny = static_cast<size_t>(std::max(r_max/Gy.norm(), 1.));
    size_t nz = static_cast<size_t>(std::max(r_max/Gz.norm(), 1.));
    this->generate_mesh(nx, ny, nz);
}

void K_mesh::generate_mesh(const size_t N)
{
    const GSL::Vector Gx = r_lattice[0], Gy = r_lattice[1], Gz = r_lattice[2];
    double vol = GSL::dot(Gx, GSL::cross(Gy, Gz));
    double Nd = static_cast<double>(N);
    size_t nx = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gx.norm()));
    size_t ny = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gy.norm()));
    size_t nz = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gz.norm()));

    this->generate_mesh(nx, ny, nz);
}