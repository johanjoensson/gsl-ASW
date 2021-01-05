#ifndef NUMERICAL_MESH_LIB_INTEGRATION_H
#define NUMERICAL_MESH_LIB_INTEGRATION_H
#include "numerical-mesh.h"
#include <functional>
#include <type_traits>

#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar steps = (end - start)/h;
    for(Scalar x = 1; x <= steps; x++){
        res +=  integrand(start + (x - 1)*h)*mesh.dr(start + (x - 1)*h) +
                integrand(start + x)*mesh.dr(start + x);
    }
    return res*h/2;
}

#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const std::function<Scalar(const Scalar&)>& integrand)
{
    return trapezoidal_integral(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), 1., integrand);
}


#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar simpson_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar steps = (end - start)/h;
    for(Scalar x = 1; x <= steps/2; x++){
        res +=  integrand(start + 2*(x - 1)*h)*mesh.dr(start + 2*(x - 1)*h)
                 + 4*integrand(start + (2*x - 1)*h)*mesh.dr(start + (2*x - 1)*h)
                 + integrand(start + 2*x*h)*mesh.dr(start + 2*x*h);
    }
    return res*h/3;
}

#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar simpson_integral(const Mesh_base<1, Scalar>& mesh, const std::function<Scalar(const Scalar&)>& integrand)
{
    return simpson_integral(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), 1., integrand);
}

#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
std::array<Scalar, K> end_point_corrections()
{
    return {9./24, 28./24, 23./24};
}

#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
Scalar corrected_trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const Scalar& h, const std::function<Scalar(const Scalar&)>& integrand)
{
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar steps = (end - start)/h;
    if(steps < 2*K){
        return 0;
    }
    auto weights = end_point_corrections<Scalar, K>();
    Scalar res = 0;
    for(Scalar i = 0; i < K; i++){
        res += weights[i]*integrand(start + i*h)*mesh.dr(start + i*h);
        res += weights[i]*integrand(end - i*h)*mesh.dr(end - i*h);
    }
    for(Scalar x = K; x <= steps - K; x++){
        res += integrand(start + x*h)*mesh.dr(start + x*h);
    }
    return res*h;
}

#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
Scalar corrected_trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const std::function<Scalar(const Scalar&)>& integrand)
{
    return corrected_trapezoidal_integral(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), 1., integrand);
}


#endif  //NUMERICAL_MESH_LIB_INTEGRATION_H
