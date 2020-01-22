#include <cmath>
#include <gsl/gsl_math.h>
#include "structure_const.h"
#include "GSLpp/special_functions.h"
#include "gaunt.h"
#include "spherical_fun.h"

Structure_constant::Structure_constant(int l_int_n, double kappa_n, lm l1_n, lm l2_n)
	: l_int(l_int_n), l1(l1_n), l2(l2_n), kappa(kappa_n)
{}

Structure_constant::Structure_constant(int l_int_n, lm l1_n, lm l2_n)
	: Structure_constant(l_int_n, std::sqrt(0.015), l1_n, l2_n)
{}

Structure_constant::Structure_constant(lm l1_n, lm l2_n)
	: Structure_constant(3, std::sqrt(0.015), l1_n, l2_n)
{}

// Structure_constant::Structure_constant()
//  : l_int(), l1(), l2(), kappa()
// {}

double Structure_constant::operator()(const GSL::Vector& r) const
{
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	double r_norm = r.norm<double>();
	int c = 1;
	if (l2.l % 2 != 0){
		c = -1;
	}

	GSL::Result sum, m_sum, a;
	Hankel_function H(lm {0, 0});
	for (int lpp = 0; lpp <= l_int; lpp++){
		H.set_l(lm{lpp, 0});
		m_sum = GSL::Result();
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l2, l1, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-15){
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r);
			}
		}

		sum += c*H(kappa*r_norm)*m_sum;
		c *= -1;
	}

	return 4*M_PI*k_fac*sum.val;
}

double Structure_constant::dot(const GSL::Vector& r) const
{
	double k_fac = GSL::pow_int(kappa, l1.l + l2.l);

	double r_norm = r.norm<double>();
	int c = 1;
	if (l2.l % 2 != 0){
		c = -1;
	}

	GSL::Result m_sum, a, tmp, sum;
	Hankel_function H(lm {0, 0});
	for (int lpp = 0; lpp <= l_int; lpp++){
		H.set_l(lm{lpp, 0});
		m_sum = GSL::Result();
		tmp = GSL::Result();
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l2, l1, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-15){
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r);
			}
		}

		tmp.val = (l1.l + l2.l - lpp)/(-kappa*kappa)*H(kappa*r_norm);
		H.set_l(lm {lpp - 1, 0});
		tmp += r_norm * H(kappa*r_norm);
		sum += c*tmp*m_sum;
		c *= -1;
	}
	return (2*M_PI*k_fac*sum).val;
}
std::ostream& operator << ( std::ostream& os, const Structure_constant& B)
{
	os << "B[(" << B.l1.l << ", " << B.l1.m << "), (" << B.l2.l << ", "
	<< B.l2.m << ")]";
	return os;
}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(
	int l_int_n, double kappa_n, const Crystal_t<3, Atom>& c_n, lm l1_n, lm l2_n)
	: l_int(l_int_n), kappa(kappa_n), c(c_n), l1(l1_n), l2(l2_n)
{}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(
	int l_int_n, const Crystal_t<3, Atom>& c_n, lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(l_int_n, std::sqrt(0.015), c_n, l1_n, l2_n)
{}

Bloch_summed_structure_constant::Bloch_summed_structure_constant(const Crystal_t<3, Atom>& c_n,
	lm l1_n, lm l2_n)
	: Bloch_summed_structure_constant(3, std::sqrt(0.015), c_n, l1_n, l2_n)
{}

GSL::Complex Bloch_summed_structure_constant::operator()(const GSL::Vector& tau,
	const GSL::Vector& kp) const
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

GSL::Complex Bloch_summed_structure_constant::dot(
	const GSL::Vector& tau, const GSL::Vector& kp) const
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
			if (std::abs(a.val) > 1E-15){
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
