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

GSL::Result cubic_harmonic(lm l, GSL::Vector& r)
{
	int m_eff = l.m;
	int c = 1;
	if (l.m < 0){
		m_eff = -l.m;
	}
	if(l.m % 2 != 0){
		c = -1;
	}
	double x = r[0];
	double y = r[1];
	double z = r[2];
	GSL::Result theta, phi, cos_theta, res;

	double r_norm = r.norm();

	theta.val = (GSL::arccos(z/r_norm)).re;
	theta.err = 0.;
	phi.val = GSL::Complex(x, y).arg();
	phi.err = 0.;
	cos_theta = GSL::cos(theta);


	if(l.m > 0){
		res = GSL::cos(m_eff*phi)*sqrt(2.);
	}else if (l.m == 0){
		res.val = 1.;
		res.err = 0;
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
