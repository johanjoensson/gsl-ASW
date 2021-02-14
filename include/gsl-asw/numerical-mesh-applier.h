#ifndef NUMERICAL_MESH_APPLIER_H
#define NUMERICAL_MESH_APPLIER_H
#include <numerical-mesh.h>

#if __cplusplus >= 202002L
template<size_t Dim, floating_point Scalar = double>
#else
template<size_t Dim, class Scalar = double>
#endif
class applier<Mesh_base<Dim, Scalar>> {
};

#endif // NUMERICAL_MESH_APPLIER_H
