#include <cmath>
#include <gsl/gsl_math.h>
#include "structure_const.h"
#include "../../GSL-lib/src/special_functions.h"
#include "gaunt.h"
#include "spherical_fun.h"

// This needs to be really thoroghly checked!
// 'cause right now it is definitely wrong!'
Structure_constant::Structure_constant(int l_low, int l_int, double kappa, lm l1, lm l2, GSL::Vector& r)
	: l_int(l_int), l_low(l_low), l1(l1), l2(l2), kappa(kappa), r(r)
{
	this->r.copy(r);

	double k_fac = GSL::pow_int(kappa, l1.l + l2.l + 1);

	double r_norm = r.norm();
	int c = 1;
	if (l1.l % 2 != 0){
		c = -1;
	}

	// Sum over all intermediate angular momenta
	// Calculate both value and energy derivative of structure constant
	GSL::Result sum, m_sum, a, tmp, d_sum;
	for (int lpp = 0; lpp <= l_int; lpp++){
		m_sum = GSL::Result();
		tmp = GSL::Result();
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l1, l2, lm {lpp, mpp});
			if (abs(a.val) > 1E-16){
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r);
			}
		}

		tmp = (l1.l + l2.l - lpp)/(kappa*kappa)*real_spherical_hankel(lm {lpp, 0}, kappa*r_norm);
		if(lpp > 0){
			tmp += r_norm/kappa * real_spherical_hankel(lm {lpp - 1, 0}, kappa*r_norm);
		}
		sum += c*real_spherical_hankel(lm {lpp, 0}, kappa*r_norm)*m_sum;
		d_sum += c*tmp*m_sum;
		c *= -1;
	}
	this->val = 4*M_PI*k_fac*sum.val;
	this->dk_val = (2*M_PI*k_fac*d_sum).val;
}

Structure_constant::Structure_constant(int l_low, int l_int, lm l1, lm l2, GSL::Vector& r)
	: Structure_constant(l_low, l_int, std::sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(int l_low, lm l1, lm l2, GSL::Vector& r)
	: Structure_constant(l_low, l_low + 1, std::sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant(lm l1, lm l2, GSL::Vector& r)
	: Structure_constant(2, 3, std::sqrt(0.015), l1, l2, r)
{
}

Structure_constant::Structure_constant()
 :l_int(), l_low(), l1(), l2(), kappa(), r(), val(), dk_val()
{

}

std::ostream& operator << ( std::ostream& os, const Structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", " << B.l2.m << ")]" << B.r <<" = " <<B.val;
	return os;
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low, int l_int, double kappa, Crystal& c, lm l1, lm l2)
	: l_low(l_low), l_int(l_int), kappa(kappa), c(c), l1(l1), l2(l2)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low, int l_int, Crystal& c, lm l1, lm l2)
	: Bloch_summed_structure_constant(l_low, l_int, std::sqrt(0.015), c, l1, l2)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low, Crystal& c, lm l1, lm l2)
	: Bloch_summed_structure_constant(l_low, l_low + 1, std::sqrt(0.015), c, l1, l2)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(Crystal& c, lm l1, lm l2)
	: Bloch_summed_structure_constant(2, 3, std::sqrt(0.015), c, l1, l2)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant()
 :l_low(), l_int(), kappa(), c(), l1(), l2()
{
}

GSL::Complex Bloch_summed_structure_constant::evaluate(const GSL::Vector& tau, const GSL::Vector& kp)
{
	Bloch_sum bs(lm {0, 0}, kappa, this->c);
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	int c = 1;
	if (l1.l % 2 != 0){
		c = -1;
	}
	GSL::Result a;
	GSL::Complex m_sum(0., 0.), sum(0., 0.);
	for (int lpp = 0; lpp <= l_int; lpp++){
		m_sum = GSL::Complex(0., 0.);
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l1, l2, lm {lpp, mpp});
			if (abs(a.val) > 1E-16){
				bs = Bloch_sum(lm {lpp, mpp}, kappa, this->c);
			    m_sum += a.val*bs.hankel_envelope(tau, kp);
			}
		}
		sum += c/(gsl_pow_int(kappa, lpp))*m_sum;
		c *= -1;
	}

	return 4*M_PI*k_fac*sum;
}
GSL::Complex Bloch_summed_structure_constant::dot_evaluate(const GSL::Vector& tau, const GSL::Vector& kp)
{

	Bloch_sum bs(lm {0, 0}, kappa, this->c);
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	int c = 1;
	if (l1.l % 2 != 0){
		c = -1;
	}
	GSL::Result a;
	GSL::Complex m_sum(0., 0.), sum(0., 0.);
	for (int lpp = 0; lpp <= l_int; lpp++){
		m_sum = GSL::Complex(0., 0.);
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l1, l2, lm {lpp, mpp});
			if (abs(a.val) > 1E-16){
				bs = Bloch_sum(lm {lpp, mpp}, kappa, this->c);
			    m_sum += a.val*(2*bs.hankel_envelope_dot(tau, kp) +
				(l1.l + l2.l - lpp)/GSL::pow_int(kappa, 2)*bs.hankel_envelope(tau, kp));
			}
		}
		sum += c/(gsl_pow_int(kappa, lpp))*m_sum;
		c *= -1;
	}
	return 2*M_PI*k_fac*sum;
}

std::ostream& operator << ( std::ostream& os, const Bloch_summed_structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", " << B.l2.m << ")]";
	return os;
}
