#include <cmath>
#include <string>
#include "spherical_fun.h"
#include "../../GSL-lib/src/complex.h"

#include <iostream>

double wronskian(Spherical_function& a, Spherical_function& b, double r)
{
    double res = 0.;
    double fl = a(r);
    a.l.l += 1;
    double flp1 = a(r);
    double gl = b(r);
    b.l.l += 1;
    double glp1 = b(r);
    res = flp1*gl - fl*glp1;

    return res;
}

double Hankel_function::operator()(const double x)
{
    GSL::Result exp = GSL::exp(-x);
    GSL::Result k;
    if(l.l >= 0){
        k = GSL::bessel_kn_scaled(l.l, x);
    }else{
        k = GSL::Result(0, 0);
    }

  return (exp*k).val;
}

double Bessel_function::operator()(const double x)
{
    GSL::Result exp = GSL::exp(-x);
    GSL::Result i;
    if(l.l >= 0){
        i = GSL::bessel_in_scaled(l.l, x);
    }else{
        i = GSL::Result(0, 0);
    }

  return (exp*i).val;
}

double Neumann_function::operator()(const double x)
{
    GSL::Result n(0,0);

    return n.val*x;
}

GSL::Result cubic_harmonic(lm l, const GSL::Vector& r)
{
	int m_eff = l.m;
	double r_norm = r.norm();
	int c = 1;
	if (l.m < 0){
		m_eff = -l.m;
	}
	if(l.m % 2 != 0){
		c = -1;
	}
	if(r_norm < 1e-14){
		GSL::Result res(0., 0.);
		return res;
	}

	GSL::Vector tmp(3);
	tmp.copy(r);
	tmp.normalize();
	double x = tmp[0];
	double y = tmp[1];
	double z = tmp[2];
	GSL::Result phi, cos_theta, res;


	phi = GSL::Result(GSL::Complex(x, y).arg(), 0.);
	cos_theta = GSL::Result(z, 0);


	if(l.m > 0){
		res = GSL::cos(m_eff*phi)*sqrt(2.);
	}else if (l.m == 0){
		res = GSL::Result(1., 0.);
	}else{
		res = GSL::sin(m_eff*phi)*sqrt(2.);
	}

	return c*GSL::legendre_sphPlm(l.l, m_eff, cos_theta.val)*res;
}

std::ostream& operator << ( std::ostream& os, const lm& l)
{
	return os << "(" << l.l << ", " << l.m << ")";
}
