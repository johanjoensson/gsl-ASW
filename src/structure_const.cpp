#include <cmath>
#include <gsl/gsl_math.h>
#include "structure_const.h"
#include "GSLpp/special_functions.h"
#include "gaunt.h"
#include "spherical_fun.h"

// This needs to be really thoroghly checked!
// 'cause right now it is almost definitely wrong!'
Structure_constant::Structure_constant(int l_low_n, int l_int_n,
	double kappa_n, lm l1_n, lm l2_n, GSL::Vector& r_n)
	: l_int(l_int_n), l_low(l_low_n), l1(l1_n), l2(l2_n), kappa(kappa_n),
	r(r_n), val(0), dk_val(0)
{
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l + 1);

	double r_norm = r.norm();
	int c = 1;
	if (l2.l % 2 != 0){
		c = -1;
	}

	// Sum over all intermediate angular momenta
	// Calculate both value and energy derivative of structure constant
	GSL::Result sum, m_sum, a, tmp, d_sum;
	Hankel_function H(lm {0, 0});
	for (int lpp = 0; lpp <= l_int; lpp++){
		H.l.l = lpp;
		m_sum = GSL::Result();
		tmp = GSL::Result();
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l2, l1, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-15){
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r);
			}
		}

		tmp.val = (l1.l + l2.l - lpp)/(-kappa*kappa)*H(kappa*r_norm);
		if(lpp > 0){
			H.l.l = lpp - 1;
			tmp += r_norm/kappa * H(kappa*r_norm);
		}
		H.l.l = lpp;
		sum += c*H(kappa*r_norm)*m_sum;
		d_sum += c*tmp*m_sum;
		c *= -1;
	}
	this->val = 4*M_PI*k_fac*sum.val;
	this->dk_val = (2*M_PI*k_fac*d_sum).val;
}

Structure_constant::Structure_constant(int l_low_n, double kappa_n, lm l1_n, lm l2_n,
	GSL::Vector& r_n)
	: Structure_constant(l_low_n, l_low_n +1, kappa_n, l1_n, l2_n, r_n)
{
}

Structure_constant::Structure_constant(int l_low_n, int l_int_n, lm l1_n, lm l2_n,
	GSL::Vector& r_n)
	: Structure_constant(l_low_n, l_int_n, std::sqrt(0.015), l1_n, l2_n, r_n)
{
}

Structure_constant::Structure_constant(int l_low_n, lm l1_n, lm l2_n, GSL::Vector& r_n)
	: Structure_constant(l_low_n, l_low_n + 1, std::sqrt(0.015), l1_n, l2_n, r_n)
{
}

Structure_constant::Structure_constant(lm l1_n, lm l2_n, GSL::Vector& r_n)
	: Structure_constant(2, 3, std::sqrt(0.015), l1_n, l2_n, r_n)
{
}

Structure_constant::Structure_constant()
 :l_int(), l_low(), l1(), l2(), kappa(), r(), val(), dk_val()
{

}

std::ostream& operator << ( std::ostream& os, const Structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", "
	<< B.l2.m << ")]" << B.r <<" = " <<B.val;
	return os;
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low_n,
	int l_int_n, double kappa_n, Crystal& c_n, lm l1_n, lm l2_n)
	: l_low(l_low_n), l_int(l_int_n), kappa(kappa_n), c(c_n), l1(l1_n), l2(l2_n)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low_n,
	double kappa_n, Crystal& c_n, lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(l_low_n, l_low_n + 1, kappa_n, c_n, l1_n, l2_n)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low_n,
	int l_int_n, Crystal& c_n, lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(l_low_n, l_int_n, std::sqrt(0.015), c_n, l1_n, l2_n)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(int l_low_n,
	Crystal& c_n, lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(l_low_n, l_low_n + 1, std::sqrt(0.015), c_n, l1_n, l2_n)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(Crystal& c_n,
	lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(2, 3, std::sqrt(0.015), c_n, l1_n, l2_n)
{
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant()
 :l_low(), l_int(), kappa(), c(), l1(), l2()
{
}

GSL::Complex Bloch_summed_structure_constant::evaluate(const GSL::Vector& tau,
	const GSL::Vector& kp)
{
	Bloch_sum bloch_sum(lm {0, 0}, kappa, this->c);
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	int sign = 1;
	GSL::Result a;
	GSL::Complex m_sum(0., 0.), sum(0., 0.);
	for (int lpp = 0; lpp <= l_int; lpp++){
		m_sum = GSL::Complex(0., 0.);
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l2, l1, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-14){
				bloch_sum = Bloch_sum(lm {lpp, mpp}, kappa, this->c);
			    m_sum += a.val*bloch_sum.hankel_envelope(tau, kp);
			}
		}
		sum += sign/(GSL::pow_int(kappa, lpp))*m_sum;
		sign = -sign;
	}

	if (l2.l % 2 == 0){
		sign = 1;
	}else{
		sign = -1;
	}
	return 4*M_PI*sign*k_fac*sum;
}

GSL::Complex Bloch_summed_structure_constant::dot_evaluate(
	const GSL::Vector& tau, const GSL::Vector& kp)
{
	Bloch_sum bloch_sum(lm {0, 0}, kappa, this->c);
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	int sign = 1;
	GSL::Result a;
	GSL::Complex m_sum(0., 0.), sum(0., 0.);
	for (int lpp = 0; lpp <= l_int; lpp++){
		m_sum = GSL::Complex(0., 0.);
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l2, l1, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-14){
				bloch_sum = Bloch_sum(lm {lpp, mpp}, kappa, this->c);
			    m_sum += a.val*(2*bloch_sum.hankel_envelope_dot(tau, kp) -
				(l1.l + l2.l - lpp)/
				GSL::pow_int(kappa, 2)*bloch_sum.hankel_envelope(tau, kp));
			}
		}
		sum += sign/(GSL::pow_int(kappa, lpp))*m_sum;
		sign = -sign;
	}
	if (l2.l % 2 == 0){
		sign = 1;
	}else{
		sign  = -1;
	}
	return 2*M_PI*sign*k_fac*sum;
}

std::ostream& operator << ( std::ostream& os,
	const Bloch_summed_structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", " << B.l2.m << ")]";
	return os;
}
