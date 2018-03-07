#include <cmath>
#include <gsl/gsl_math.h>
#include "structure_const.h"
#include "gaunt.h"
#include "spherical_fun.h"

Structure_constant::Structure_constant(int l_low, int l_int, double kappa, lm l1, lm l2, gsl_vector r)
	: l1(l1), l2(l2), r(r)
{
	this->l_low = l_low;
	this->l_int = l_int;
	this->kappa = kappa;

	double k_fac = gsl_pow_int(kappa, l1.l + l2.l + 1);
	
	double x = gsl_vector_get(&r, 0);
	double y = gsl_vector_get(&r, 1);
	double z = gsl_vector_get(&r, 2);
	double r_norm = std::sqrt(x*x + y*y + z*z);



	int c = 1;
	if (l1.l % 2 != 0){
		c = -1;
	}

	// Sum over all intermediate angular momenta
	double sum = 0., m_sum = 0., a = 0.;
	for (int lpp = 0; lpp <= l_int; lpp++){
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l1, l2, lm {lpp, mpp});
			if (a != 0){
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r);
			}
		}
		sum += c*real_spherical_hankel(lm {lpp, 0}, kappa*r_norm)*m_sum;
		c *= -1;
	}
	this->val = 4*M_PI*k_fac*sum;
}

Structure_constant::Structure_constant(int l_low, int l_int, lm l1, lm l2, gsl_vector r)
	: Structure_constant(l_low, l_int, std::sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(int l_low, lm l1, lm l2, gsl_vector r)
	: Structure_constant(l_low, l_low + 1, std::sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(lm l1, lm l2, gsl_vector r)
	: Structure_constant(2, 3, std::sqrt(0.015), l1, l2, r)
{
}


std::ostream& operator << ( std::ostream& os, const Structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", " << B.l2.m << ")](" << gsl_vector_get(&B.r,0) << ", " << gsl_vector_get(&B.r,1) << ", " << gsl_vector_get(&B.r,2) << ") = " <<B.val;
	return os;
}

