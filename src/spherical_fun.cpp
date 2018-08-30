#include <cmath>
#include <string>
#include "spherical_fun.h"
#include "../../GSL-lib/src/complex.h"

#include <iostream>

unsigned long int factorial(int n)
{
	unsigned long int res = 1;
	for (int i = 1; i <= n; i++){
		res *= i;
	}

	return res;

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
	if(r_norm < 1e-16){
		GSL::Result res(0., 0.);
		return res;
	}

	GSL::Vector tmp(3);
	tmp.copy(r);
	tmp.normalize();
	double x = tmp[0];
	double y = tmp[1];
	double z = tmp[2];
	GSL::Result theta, phi, cos_theta, res;


	theta = GSL::Result((GSL::arccos(z/r_norm)).re, 0.);
	phi = GSL::Result(GSL::Complex(x, y).arg(), 0.);
	cos_theta = GSL::cos(theta);


	if(l.m > 0){
		res = GSL::cos(m_eff*phi)*sqrt(2.);
	}else if (l.m == 0){
		res = GSL::Result(1., 0.);
	}else{
		res = GSL::sin(m_eff*phi)*sqrt(2.);
	}

	return c*GSL::legendre_sphPlm(l.l, m_eff, cos_theta.val)*res;
}

GSL::Result real_spherical_hankel(lm l, double x)
{
  GSL::Result exp = GSL::exp(-x);
  GSL::Result k = 2./M_PI * GSL::bessel_kn_scaled(l.l, x);

  return exp*k;
}

std::ostream& operator << ( std::ostream& os, const lm& l)
{
	return os << "(" << l.l << ", " << l.m << ")";
}
