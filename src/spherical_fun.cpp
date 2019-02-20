#include <cmath>
#include <string>
#include "spherical_fun.h"
#include "GSLpp/complex.h"

#include <iostream>

double wronskian(Spherical_function& a, Spherical_function& b, double r)
{
    double res = 0.;
    double fl = a(r);
    a.set_l(lm {a.l().l + 1, a.l().m});
    double flp1 = a(r);
    double gl = b(r);
    b.set_l(lm {b.l().l + 1, b.l().m});
    double glp1 = b(r);
    res = flp1*gl - fl*glp1;

    return res;
}

double Hankel_function::operator()(const double x) const
{
    GSL::Result exp = GSL::exp(-x);
    GSL::Result k;
    k = GSL::bessel_kn_scaled(l_m.l, x);
    return 2./M_PI*(exp*k).val;
}

double Bessel_function::operator()(const double x) const
{
    GSL::Result exp = GSL::exp(std::abs(x));
    GSL::Result i = GSL::bessel_in_scaled(l_m.l, x);

  return (exp*i).val;
}

double Neumann_function::operator()(const double x) const
{
    GSL::Result n(0,0);

    return n.val*x;
}

double Integral_Hankel_function::operator()(const double x) const
{
	return GSL::pow_int(2*x, l_m.l)*( I.ewald_int(l_m, x) +
	  2./M_SQRTPI *I.comp_ewald_int(l_m, x));
}

void Integral_Hankel_function::set_ewald_param(const double eta)
{
	I.set_ewald_param(eta);
}

GSL::Result cubic_harmonic(lm l, const GSL::Vector& r)
{
	int m_eff = std::abs(l.m);
    int sign = 1;
    if(m_eff % 2 == 1){
        sign = -1;
    }
	double r_norm = r.norm<double>();
	if(r_norm < 1e-14){
		GSL::Result res(0., 0.);
		return res;
	}

	GSL::Vector tmp(3);
	tmp.copy(r);
	tmp.normalize<double>();
	double x = tmp[0];
	double y = tmp[1];
	double z = tmp[2];
	GSL::Result phi, cos_theta, res;


	phi = GSL::Result(GSL::Complex(x, y).arg(), 0.);
	cos_theta = GSL::Result(z, 0);


	if(l.m >= 0){
		res = GSL::cos(m_eff*phi);
	}else if (l.m < 0){
		res = GSL::sin(m_eff*phi);
	}
    if(l.m != 0){
        res *= std::sqrt(2.);
    }
	return sign*GSL::legendre_sphPlm(l.l, m_eff, cos_theta.val)*res;
}
