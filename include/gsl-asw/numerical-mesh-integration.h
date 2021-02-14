#ifndef NUMERICAL_MESH_LIB_INTEGRATION_H
#define NUMERICAL_MESH_LIB_INTEGRATION_H

#if __cplusplus < 201103L // Require at least c++11
#error The numerical-mesh-integration library requires at least c++11!
#endif

#include <numerical-mesh.h>
#include <functional>
#include <type_traits>

/*******************************************************************************
* The __cplusplus macro is used to determine which c++ standard is used for the*
* compilation. This is "needed" because I want to use this library in projects *
* using c++11 as well as projects using c++20. This means parts of the code had*
* to be rewritten in order to enable the cool c++20 features only if c++20 is  *
* used. As a side effect this also made it possible to write code for c++17 as *
* well.                                                                        *
* Features requiring c++17:                                                    *
*  -- structured bindings i.e., auto [a,b,c,...] = ...                         *
* Features requiring c++20:                                                    *
*  -- Concepts e.g., integral and floating_point                               *
*  -- ranges-versions of algorithms                                            *
*  -- std::transform with 2 input containers                                   *
*******************************************************************************/

/*!*****************************************************************************
 * \defgroup NumMeshInt Numerical integration
 * @{
 ******************************************************************************/

//! Trapezoidal method numerical integral over numerical mesh
/*!*****************************************************************************
* Trapezoidal integration of the integrand over the mesh, from first meshpoint
* to last.
* \f[
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r \approx \frac{1}{2} \sum_{i = x_0}^{x^1}
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_x f(r(x))
* \f]
* @param mesh Numerical mesh to use for the integration.
* @param integrand Function to integrate.
* \return The approximate value of the integral.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar trapezoidal_integral(const Mesh_base<1, Scalar>& mesh,
    const std::function<Scalar(const Scalar&)>& integrand)
    noexcept
{
    return trapezoidal_integral<Scalar>(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), integrand, 1);
}

//! Compound Simpson's rule method numerical integral over numerical mesh
/*!*****************************************************************************
* Compound Simpson rule integration of the integrand over the mesh, from first
* meshpoint to last.
* \f[
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r \approx \frac{1}{3} \sum_{i = 1}^{n/2}
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + 2(x_i-1)step}f(2(x_i - 1)step)
* + 4\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + (2x_i-1)step}f((2x_i - 1)step)
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + 2x_i step}f(2x_i step)
* \f]
* @param mesh Numerical mesh to use for the integration.
* @param integrand Function to integrate.
* \return The approximate value of the integral.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar simpson_integral(const Mesh_base<1, Scalar>& mesh,
    const std::function<Scalar(const Scalar&)>& integrand)
    noexcept
{
    return simpson_integral<Scalar>(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), integrand, 1);
}

//! Trapezoidal rule, with modified end points, method numerical integral over numerical mesh
/*!*****************************************************************************
* Trapezoidal integration, with modified end point weight, of the integrand over
* the mesh, from meshpoint x0 to mesh point x1, taking steps of size step.
* \f{eqnarray*}{
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r &\approx&
* \sum_{i = 0}^{k-1} w(i)*\left(
*       f(r(x_0 + i*step))\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + i*step}
*     + f(r(x_1 - i*step))\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_1 - i*step}
* \right)\\
* &+& \sum_{x_i = x_0 + i*k}^{x1 - i*k}
*   \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_i} f(r(x_i))
* \f}
* @param mesh Numerical mesh to use for the integration.
* @param integrand Function to evaluate.
* \return The approximate value of the integrand.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
Scalar corrected_trapezoidal_integral(const Mesh_base<1, Scalar>& mesh,
    const std::function<Scalar(const Scalar&)>& integrand)
    noexcept
{
    return corrected_trapezoidal_integral<Scalar, K>(mesh, 0., static_cast<Scalar>(mesh.dim() - 1), integrand, 1);
}

namespace{
#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
std::array<Scalar, K> end_point_corrections()
    noexcept
{
    return {9./24, 28./24, 23./24};
}
}

//! Trapezoidal method numerical integral over part of numerical mesh
/*!*****************************************************************************
* Trapezoidal integration of the integrand over the mesh, from meshpoint x0
* to mesh point x1, taking steps of size step.
* \f[
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r \approx \frac{1}{2} \sum_{i = x_0}^{x^1}
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_i} f(r(x_i)), \, x_i = r(x_0 + i*step)
* \f]
* @param mesh Numerical mesh to use for the integration.
* @param x0 Lower end of the integral (in mesh point indices).
* @param x1 Upper end of the integral (in mesh point indices).
* @param step Step size (in mesh point units).
* @param integrand Function to evaluate.
* \return The approximate value of the integrand.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const std::function<Scalar(const Scalar&)>& integrand, const Scalar& step = 1)
    noexcept
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);

    for(Scalar x = start + step; x <= end; x += step){
        res +=  integrand(x - step)*mesh.dr(x - step) +
                integrand(x)*mesh.dr(x);
    }

    return res*step/2;
}

//! Compound Simpson's rule method numerical integral over part of numerical mesh
/*!*****************************************************************************
* Composite Simpsons rule integration of the integrand over the mesh, from
* meshpoint x0 to mesh point x1, taking steps of size step.
* \f[
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r \approx \frac{1}{3} \sum_{i = 1}^{n/2}
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + 2(x_i-1)step}f(2(x_i - 1)step)
* + 4\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + (2x_i-1)step}f((2x_i - 1)step)
*  \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + 2x_i step}f(2x_i step)
* \f]
* @param mesh Numerical mesh to use for the integration.
* @param x0 Lower end of the integral (in mesh point indices).
* @param x1 Upper end of the integral (in mesh point indices).
* @param step Step size (in mesh point units).
* @param integrand Function to evaluate.
* \return The approximate value of the integrand.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar>
#else
template<class Scalar>
#endif
Scalar simpson_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const std::function<Scalar(const Scalar&)>& integrand, const Scalar& step = 1)
    noexcept
{
    Scalar res = 0;
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    for(Scalar x = start + step; x <= (end - start)/2; x += step){
        res +=  integrand(2*x - 2*step)*mesh.dr(2*x - 2*step)
                 + 4*integrand(2*x - step)*mesh.dr(2*x - step)
                 + integrand(2*x)*mesh.dr(2*x);
    }
    return res*step/3;
}

//! Trapezoidal rule, with modified end points, method numerical integral over part of numerical mesh
/*!*****************************************************************************
* Trapezoidal integration, with modified end point weight, of the integrand over
* the mesh, from meshpoint x0 to mesh point x1, taking steps of size step.
* \f{eqnarray*}{
*  \int_{r(x_0)}^{r(x^1)}\, f(r(x)) \mathrm{d}r &\approx&
* \sum_{i = 0}^{k-1} w(i)*\left(
*       f(r(x_0 + i*step))\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_0 + i*step}
*     + f(r(x_1 - i*step))\left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_1 - i*step}
* \right)\\
* &+& \sum_{x_i = x_0 + i*k}^{x1 - i*k}
*   \left.\frac{\mathrm{d}r}{\mathrm{d}x}\right|_{x = x_i} f(r(x_i))
* \f}
* @param mesh Numerical mesh to use for the integration.
* @param x0 Lower end of the integral (in mesh point indices).
* @param x1 Upper end of the integral (in mesh point indices).
* @param step Step size (in mesh point units).
* @param integrand Function to evaluate.
* \return The approximate value of the integrand.
*******************************************************************************/
#if __cplusplus >= 202002L
template<std::floating_point Scalar, size_t K = 3>
#else
template<class Scalar, size_t K = 3>
#endif
Scalar corrected_trapezoidal_integral(const Mesh_base<1, Scalar>& mesh, const Scalar& x0,
    const Scalar& x1, const std::function<Scalar(const Scalar&)>& integrand, const Scalar& step = 1)
    noexcept
{
    Scalar start = std::min(x0, x1), end = std::max(x0, x1);
    Scalar steps = (end - start)/step;
    if(steps < 2*K){
        return 0;
    }
    auto weights = end_point_corrections<Scalar, K>();
    Scalar res = 0;
    for(Scalar i = 0; i < K; i++){
        res += weights[i]*integrand(start + i*step)*mesh.dr(start + i*step);
        res += weights[i]*integrand(end - i*step)*mesh.dr(end - i*step);
    }
    for(Scalar x = start + K*step; x <= end - K*step; x += step){
        res += integrand(x)*mesh.dr(x);
    }
    return res*step;
}

/*!*****************************************************************************
*@}
*******************************************************************************/

#endif  //NUMERICAL_MESH_LIB_INTEGRATION_H
