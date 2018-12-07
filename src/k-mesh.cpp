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
    for(size_t n1 = 0; n1 < nx; n1++){
        for(size_t n2 = 0; n2 < ny; n2++){
            for(size_t n3 = 0; n3 < nz; n3++){
                tmp = (2.*n1 - nx - 1)/(2.*nx)*Gx;
                tmp += (2.*n2 - ny - 1)/(2.*ny)*Gy;
                tmp += (2.*n3 - nz - 1)/(2.*nz)*Gz;
                k_points.push_back(GSL::Vector(tmp));
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
    size_t nx = static_cast<size_t>(std::round(std::cbrt(N/vol)*Gx.norm()));
    size_t ny = static_cast<size_t>(std::round(std::cbrt(N/vol)*Gy.norm()));
    size_t nz = static_cast<size_t>(std::round(std::cbrt(N/vol)*Gz.norm()));

    this->generate_mesh(nx, ny, nz);
}
