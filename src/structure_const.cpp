#include <cmath>
#include <gsl/gsl_math.h>
#include <algorithm>
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
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r/r.norm<double>());
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
			    m_sum += a*cubic_harmonic(lm {lpp, mpp}, r/r);
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

std::unordered_map<Bloch_summed_structure_constant::Key, GSL::Complex> Bloch_summed_structure_constant::values_m{};
std::unordered_map<Bloch_summed_structure_constant::Key, GSL::Complex> Bloch_summed_structure_constant::dot_values_m{};

GSL::Complex Bloch_summed_structure_constant::operator()(
	const Crystal_t<3, Atom>& c, const lm lint, const lm& l, const lm& lp, const double kappa,
	const GSL::Vector& tau,	const GSL::Vector& kp) const
{
	auto it = values_m.find({l, lp, kappa, tau, kp});
	if(it != values_m.end()){
		return it->second;
	}else{
		double k_fac = GSL::pow_int(kappa, l.l + lp.l);

		int sign = 1;
		GSL::Result a;
		GSL::Complex m_sum(0., 0.), sum(0., 0.);
		Bloch_sum D;
		for (int lpp = 0; lpp <= lint.l; lpp++){
			m_sum = GSL::Complex(0., 0.);
			for (int mpp = -lpp; mpp <= lpp; mpp++){
				a = gaunt(lp, l, lm {lpp, mpp});
				if (std::abs(a.val) > 1E-14){
				    m_sum += a.val*D(lm {lpp, mpp}, kappa, c, tau, kp);
				}
			}
			sum += GSL::pow_int(-1/kappa, lpp)*m_sum;
		}

		if (lp.l % 2 == 0){
			sign = 1;
		}else{
			sign = -1;
		}
		values_m.insert({{l, lp, kappa, tau, kp},{4*M_PI*sign*k_fac*sum}});
		return 4*M_PI*sign*k_fac*sum;
	}
}

GSL::Complex Bloch_summed_structure_constant::dot(
	const Crystal_t<3, Atom>& c, const lm lint, const lm& l, const lm& lp, const double kappa,
	const GSL::Vector& tau, const GSL::Vector& kp) const
{
	auto it = dot_values_m.find({l, lp, kappa, tau, kp});
	if(it != dot_values_m.end()){
		return it->second;
	}else{
		double k_fac = GSL::pow_int(kappa, l.l + lp.l);

		int sign = 1;
		GSL::Result a;
		GSL::Complex m_sum(0., 0.), sum(0., 0.);
		Bloch_sum D;
		for (int lpp = 0; lpp <= lint.l; lpp++){
			m_sum = GSL::Complex(0., 0.);
			for (int mpp = -lpp; mpp <= lpp; mpp++){
				a = gaunt(lp, l, lm {lpp, mpp});
				if (std::abs(a.val) > 1E-15){
				    m_sum += a.val*
					(2*D.dot(lm {lpp, mpp}, kappa, c, tau, kp) -
					(l.l + lp.l - lpp)/
					GSL::pow_int(kappa, 2)*D(lm {lpp, mpp}, kappa, c, tau, kp));
				}
			}
			sum += sign/(GSL::pow_int(kappa, lpp))*m_sum;
			sign = -sign;
		}
		if (lp.l % 2 == 0){
			sign = 1;
		}else{
			sign  = -1;
		}
		dot_values_m.insert({{l, lp, kappa, tau, kp},{2*M_PI*sign*k_fac*sum}});
		return 2*M_PI*sign*k_fac*sum;
	}
}
