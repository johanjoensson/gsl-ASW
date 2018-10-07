#include "envelope_fun.h"
#include "spherical_fun.h"

Envelope_function::Envelope_function()
 : center(), l(), kappa()
{}

Envelope_function::Envelope_function(const GSL::Vector center, const lm l, const double kappa)
 : center(center), l(l), kappa(kappa)
{}

double Envelope_function::operator()(GSL::Vector r)
{
    double x = (r - center).norm();
    return (this->barred_fun(x)*cubic_harmonic(l,r)).val;
}

double Envelope_Hankel::barred_fun(double x)
{
    Hankel_function h(l);
    return GSL::pow_int(kappa, l.l+1)*h(kappa*x);
}

double Envelope_Bessel::barred_fun(double x)
{
    Bessel_function j(l);
    return GSL::pow_int(kappa, -l.l)*j(kappa*x);
}

double Envelope_Neumann::barred_fun(double x)
{
    Neumann_function n(l);
    return GSL::pow_int(kappa, l.l+1)*n(kappa*x);
}
