#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "spherical_fun.h"

unsigned long int factorial(int n)
{
	unsigned long int res = 1;
	for (int i = 1; i <= n; i++){
		res *= i;
	}

	return res;

}

double cubic_harmonic(lm l, gsl_vector r)
{
	int m_eff = l.m;
	int c = 1;
	if (l.m % 2 != 0){
		c = -1;
	}
	if (l.m < 0){
		m_eff = -l.m;
	}
	double fac = double(factorial(l.l - m_eff))/factorial(l.l - m_eff);
	double x = gsl_vector_get(&r, 0);
	double y = gsl_vector_get(&r, 1);
	double z = gsl_vector_get(&r, 2);
	double theta, phi, cos_theta, cos_phi, sin_phi, res = 0;

	double r_norm = std::sqrt(x*x + y*y + z*z);

	theta = gsl_complex_arccos_real(z/r_norm).dat[0];
	phi = gsl_complex_arctan(gsl_complex_rect(y/x, 0)).dat[0];

	cos_theta = gsl_sf_cos(theta);

	if(l.m > 0){
		cos_phi = gsl_sf_cos(m_eff*phi);
		res = c*std::sqrt(2)*fac*gsl_sf_legendre_sphPlm(l.l, m_eff, cos_theta)*cos_phi;
	}else if (l.m == 0){
		cos_phi = gsl_sf_cos(m_eff*phi);
		res = fac*gsl_sf_legendre_sphPlm(l.l, m_eff, cos_theta)*cos_phi;
	}else{
		sin_phi = gsl_sf_sin(m_eff*phi);
		res = c*std::sqrt(2)*fac*gsl_sf_legendre_sphPlm(l.l, m_eff, cos_theta)*sin_phi;;
	}
		
	return res;
}

double real_spherical_hankel(lm l, double x)
{
  double exp = gsl_sf_exp(-x);
  double k = gsl_sf_bessel_kl_scaled(l.l, x);

  return exp*k;
}

std::ostream& operator << ( std::ostream& os, const lm& l) 
{
	return os << "(" << l.l << ", " << l.m << ")";
}
