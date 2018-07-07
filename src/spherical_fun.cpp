#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <string>
#include "spherical_fun.h"
#include "../../GSL-lib/src/complex.h"
#include "../../GSL-lib/src/special_functions.h"


unsigned long int factorial(int n)
{
	unsigned long int res = 1;
	for (int i = 1; i <= n; i++){
		res *= i;
	}

	return res;

}

double cubic_harmonic(lm l, GSL::Vector& r)
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
	double theta = 0, phi = 0, cos_theta = 0, res = 0;

	double r_norm = r.norm();

	theta = (GSL::arccos(z/r_norm)).re;
	phi = GSL::Complex(x, y).arg();
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
  double k = 2./M_PI * std::get<0>(GSL::bessel_kn_scaled(l.l, x));

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
