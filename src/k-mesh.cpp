#include <k-mesh.h>
#include <iostream>
#include <cmath>

K_mesh::K_mesh()
 : r_lattice(), k_points()
{}

K_mesh::K_mesh(const GSL::Matrix& r_lat)
 : r_lattice(r_lat), k_points()
{}

/*******************************************************************************
* Monkhorst-Pack scheme for generatign K-mesh                                  *
*******************************************************************************/
void K_mesh::generate_mesh(const size_t nx, const size_t ny, const size_t nz)
{
    const GSL::Vector Gx(r_lattice[0]), Gy(r_lattice[1]), Gz(r_lattice[2]);
    GSL::Vector kp, tmp;
    double nxd = static_cast<double>(nx);
    double nyd = static_cast<double>(ny);
    double nzd = static_cast<double>(nz);

    for(size_t n1 = 0; n1 < nx; n1++){
        for(size_t n2 = 0; n2 < ny; n2++){
            for(size_t n3 = 0; n3 < nz; n3++){
                k_points.push_back(
                    r_lattice[0]*(static_cast<double>(2*n1) - nxd - 1)/(static_cast<double>(2*nx))
                +   r_lattice[1]*(static_cast<double>(2*n2) - nyd - 1)/(static_cast<double>(2*ny))
                +   r_lattice[2]*(static_cast<double>(2*n3) - nzd - 1)/(static_cast<double>(2*nz))
                );
            }
        }
    }
}

void K_mesh::generate_mesh(const double r_max)
{
    const GSL::Vector Gx = r_lattice[0], Gy = r_lattice[1], Gz = r_lattice[2];
    size_t nx = static_cast<size_t>(std::max(r_max/Gx.norm<double>(), 1.));
    size_t ny = static_cast<size_t>(std::max(r_max/Gy.norm<double>(), 1.));
    size_t nz = static_cast<size_t>(std::max(r_max/Gz.norm<double>(), 1.));
    this->generate_mesh(nx, ny, nz);
}

void K_mesh::generate_mesh(const size_t N)
{
    const GSL::Vector Gx = r_lattice[0], Gy = r_lattice[1], Gz = r_lattice[2];
    double vol = Gx.dot(Gy.cross(Gz));
    double Nd = static_cast<double>(N);
    size_t nx = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gx.norm<double>()));
    size_t ny = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gy.norm<double>()));
    size_t nz = static_cast<size_t>(std::round(std::cbrt(Nd/vol)*Gz.norm<double>()));

    this->generate_mesh(nx, ny, nz);
}

void K_mesh::generate_mesh(const std::vector<GSL::Vector>& path, const size_t N_steps)
{
    double dn = 1./static_cast<double>(N_steps);
	for(size_t i = 0; i < path.size() - 1; i++){
		for(unsigned int n = 0; n < N_steps; n++ ){
			k_points.push_back((path[i] -
					n*dn*(path[i] - path[i + 1]))*r_lattice);
		}
	}
}
