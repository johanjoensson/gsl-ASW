#ifndef K_MESH_H
#define K_MESH_H

#include <vector>
#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"

class K_mesh {
private:
    GSL::Matrix r_lattice;
public:
    std::vector<GSL::Vector> k_points;

    K_mesh() = default;
    K_mesh(GSL::Matrix&& r_lat);
    void generate_mesh(const size_t nx, const size_t ny, const size_t nz);
    void generate_mesh(const double r_max);
    void generate_mesh(const size_t N);
    void generate_mesh(const std::vector<GSL::Vector>& path, const size_t N_steps);

};
#endif
