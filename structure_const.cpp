#include <math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include "structure_const.h"
#include "gaunt.h"
#include "spherical_fun.h"

Structure_constant::Structure_constant(int l_low, int l_int, double kappa, lm l1, lm l2, gsl_vector r)
{
	this->l_low = l_low;
	this->l_int = l_int;
	this->kappa = kappa;
	this->l1 = l1;
	this->l2 = l2;
	this->r = r;

	double k_fac = gsl_pow_int(kappa, l1.l + l2.l + 1);
	
	double x = gsl_vector_get(&r, 0);
	double y = gsl_vector_get(&r, 1);
	double z = gsl_vector_get(&r, 2);
	double r_norm = sqrt(x*x + y*y + z*z);



	int c = 1;
	if (l1.l % 2 != 0){
		c = -1;
	}

	double sum, m_sum, a = 0;
	for (int lpp = 0; lpp <= l_int; lpp++){
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = real_gaunt(l1.l, l2.l, lpp, l1.m, l2.m, mpp);
			if (a != 0){
			    m_sum += a*cubic_harmonic(lpp, mpp, r);
			}
		}
		sum += c*real_spherical_hankel(lpp, kappa*r_norm)*m_sum;
		c *= -1;
	}
	this->val = 4*M_PI*k_fac*sum;
}

Structure_constant::Structure_constant(int l_low, int l_int, lm l1, lm l2, gsl_vector r)
	: Structure_constant(l_low, l_int, sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(int l_low, lm l1, lm l2, gsl_vector r)
	: Structure_constant(l_low, l_low + 1, sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(lm l1, lm l2, gsl_vector r)
	: Structure_constant(2, 3, sqrt(0.015), l1, l2, r)
{
}
