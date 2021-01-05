#include <cmath>
#include <string>
#include <spherical_fun.h>
#include <GSLpp/complex.h>

#include <iostream>

double wronskian(Spherical_function a, Spherical_function b, double r)
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
    return GSL::bessel_yn(l_m.l, x).val;
}

double Integral_Hankel_function::operator()(const double kappa, const double eta, const double x) const
{
	return GSL::pow_int(2*x, l_m.l)*( I.ewald_int(kappa, eta, l_m, x) +
	  2./M_SQRTPI*I.comp_ewald_int(kappa, eta, l_m, x));
}

GSL::Result cubic_harmonic(const lm& l, GSL::Vector&& r)
{
    if(l.l == 0 && l.m == 0){
        return GSL::Result(1./std::sqrt(4*M_PI), 0);
    }
    if(std::abs(r.norm()) < 1e-15){
        return GSL::Result(0, 0);
    }

    return GSL::pow_int(r.norm(), l.l)*cubic_harmonic(l.l, l.m,
                            r[2]/r.norm(), GSL::Complex(r[0], r[1]).arg());
}

GSL::Result cubic_harmonic(const lm& l, const GSL::Vector& r)
{
    if(l.l == 0){
        return GSL::Result(1./std::sqrt(4*M_PI), 0);
    }

    if(std::abs(r.norm()) < 1e-15){
        return GSL::Result(0, 0);
    }
	return GSL::pow_int(r.norm(), l.l)*cubic_harmonic(l.l, l.m,
                            r[2]/r.norm(), GSL::Complex(r[0], r[1]).arg());
}

GSL::Result cubic_harmonic(const int l, const int m, const double cos_theta, const double phi)
{
	int m_eff = std::abs(m);
    int sign = 1;
    if(m_eff % 2 == 1){
        sign = -1;
    }
    double fac = 1;
	if(m > 0){
		fac = GSL::cos(m_eff*phi).val;
	}else if (m < 0){
		fac = GSL::sin(m_eff*phi).val;
	}
    if(m != 0){
        fac *= std::sqrt(2.);
    }
	return sign*GSL::legendre_sphPlm(l, m_eff, cos_theta)*fac;
}
