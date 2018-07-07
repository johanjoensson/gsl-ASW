#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <string>
#include "spherical_fun.h"

//#include <iostream>

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
	if (l.m < 0){
		m_eff = -l.m;
	}
	if(l.m % 2 != 0){
		c = -1;
	}
	double x = gsl_vector_get(&r, 0);
	double y = gsl_vector_get(&r, 1);
	double z = gsl_vector_get(&r, 2);
	double theta = 0, phi = 0, cos_theta = 0, res = 0;

	double r_norm = std::sqrt(x*x + y*y + z*z);

	theta = gsl_complex_arccos_real(z/r_norm).dat[0];
	phi = gsl_complex_arg(gsl_complex_rect(x, y));
	cos_theta = gsl_sf_cos(theta);

	if(l.m > 0){
		res = gsl_sf_cos(m_eff*phi)*sqrt(2);
	}else if (l.m == 0){
		res = 1.;
	}else{
		res = gsl_sf_sin(m_eff*phi)*sqrt(2);
	}
	return c*gsl_sf_legendre_sphPlm(l.l, m_eff, cos_theta)*res;
}

double real_spherical_hankel(lm l, double x)
{
  double exp = exp_gsl(-x);
  double k = 2./M_PI * gsl_sf_bessel_kl_scaled(l.l, x);

  return exp*k;
}

double exp_gsl(double x)
{
	gsl_sf_result res;
	res.val = 0.;
	res.err = 0.;
	int stat = 0;
	stat = gsl_sf_exp_e(x, &res);
	if(stat == GSL_EUNDRFLW){
		res.val = 0.;
	}else if(stat){
		std::string error_str =   gsl_strerror(stat);
		throw std::runtime_error("Error in exponential.\nGSL error: "
		+ error_str);
	}

	return res.val;
}

std::ostream& operator << ( std::ostream& os, const lm& l)
{
	return os << "(" << l.l << ", " << l.m << ")";
}

gsl_complex operator + (gsl_complex a, gsl_complex b)
{
	return gsl_complex_add(a, b);
}

std::ostream& operator << ( std::ostream& os, const gsl_complex& a)
{
	if(GSL_IMAG(a) < 0){
		return os << GSL_REAL(a) << " " << GSL_IMAG(a) << "i";
	}else{
		return os << GSL_REAL(a) << "  + " << GSL_IMAG(a) << "i";
	}
}
gsl_complex operator - (gsl_complex a, gsl_complex b)
{
	return gsl_complex_sub(a, b);
}
gsl_complex operator * (gsl_complex a, gsl_complex b)
{
	return gsl_complex_mul(a, b);
}

gsl_complex operator / (gsl_complex a, gsl_complex b)
{
	return gsl_complex_div(a, b);
}

gsl_complex operator + (gsl_complex a, double b)
{
	return gsl_complex_add_real(a, b);
}
gsl_complex operator - (gsl_complex a, double b)
{
	return gsl_complex_sub_real(a, b);
}
gsl_complex operator * (gsl_complex a, double b)
{
	return gsl_complex_mul_real(a, b);
}
gsl_complex operator / (gsl_complex a, double b)
{
	return gsl_complex_div_real(a, b);
}

gsl_complex operator + (double a, gsl_complex b)
{
	return gsl_complex_add_real(b, a);
}
gsl_complex operator - (double a, gsl_complex b)
{
	return gsl_complex_sub_real(b, a);
}
gsl_complex operator * (double a, gsl_complex b)
{
	return gsl_complex_mul_real(b, a);
}
gsl_complex operator / (double a, gsl_complex b)
{
	return gsl_complex_mul_real(gsl_complex_inverse(b), a);
}

gsl_vector& operator + (gsl_vector& a, gsl_vector& b)
{
	gsl_vector *res = gsl_vector_alloc(a.size);
	gsl_vector_memcpy(res, &a);
	gsl_vector_add(res, &b);
	return *res;
}
