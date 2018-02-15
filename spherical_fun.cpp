#include <math.h>
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

double cubic_harmonic(int l, int m, gsl_vector r)
{
	int m_eff = m;
	int c = 1;
	if (m % 2 != 0){
		c = -1;
	}
	double fac = 1;
	if (m < 0){
		fac = c*double(factorial(l+m))/factorial(l-m);
		m_eff = -m;
	}
	double x = gsl_vector_get(&r, 0);
	double y = gsl_vector_get(&r, 1);
	double z = gsl_vector_get(&r, 2);
	double theta, phi, cos_theta, cos_phi = 0;

	double r_norm = sqrt(x*x + y*y + z*z);

	theta = gsl_complex_arccos_real(z/r_norm).dat[0];
	phi = gsl_complex_arctan(gsl_complex_rect(y/x, 0)).dat[0];

	cos_theta = gsl_sf_cos(theta);
	cos_phi = gsl_sf_cos(m_eff*phi);

	return c*fac*gsl_sf_legendre_sphPlm(l, m_eff, cos_theta)*cos_phi;
}

double real_spherical_hankel(int l, double x)
{
  double exp = gsl_sf_exp(-x);
  double k = gsl_sf_bessel_kl_scaled(l, x);

  return exp*k;
}
