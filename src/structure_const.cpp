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

void Bloch_summed_structure_constant::Container::add(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau,	const GSL::Vector& kp)
{
	Bloch_summed_structure_constant B;
	values_m.insert({{l, lp, kappa, tau, kp},B.calc(c, lint, lp, l, kappa, tau, kp, D_m)});
	dot_values_m.insert({{l, lp, kappa, tau, kp},B.dot(c, lint, lp, l, kappa, tau, kp, D_m)});
}

GSL::Complex Bloch_summed_structure_constant::Container::get(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau,	const GSL::Vector& kp)
{
/*
	// Check if value has already been calculated
	auto it = values_m.find({lp, l, kappa, tau, kp});
	if(it != values_m.end()){
		return it->second;
	}
	it = values_m.find({l, lp, kappa, tau, kp});
	if(it != values_m.end()){
		return (l.l - lp.l) % 2 == 0 ? it->second : -it->second;
	}
	it = values_m.find({lp, l, kappa, -tau, -kp});
	if(it != values_m.end()){
		return (l.l - lp.l) % 2 == 0 ? it->second : -it->second;
	}
	it = values_m.find({lp, l, kappa, tau, -kp});
	if(it != values_m.end()){
		return it->second.conjugate();
	}
	// Value not already calculated, calculate!
	this->add(c, lint, lp, l, kappa, tau, kp);
*/
	const Bloch_summed_structure_constant B;
	return B.calc(c, lint, lp, l, kappa, tau, kp, D_m);
}

GSL::Complex Bloch_summed_structure_constant::Container::get_dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa,
	const GSL::Vector& tau,	const GSL::Vector& kp)
{
/*
	// Check if value has already been calculated
	auto it = dot_values_m.find({lp, l, kappa, tau, kp});
	if(it != dot_values_m.end()){
		return it->second;
	}
	it = dot_values_m.find({l, lp, kappa, tau, kp});
	if(it != dot_values_m.end()){
		return (l.l - lp.l) % 2 == 0 ? it->second : -it->second;
	}
	it = dot_values_m.find({lp, l, kappa, -tau, -kp});
	if(it != dot_values_m.end()){
		return (l.l - lp.l) % 2 == 0 ? it->second : -it->second;
	}
	it = dot_values_m.find({lp, l, kappa, tau, -kp});
	if(it != dot_values_m.end()){
		return it->second.conjugate();
	}
	// Value not already calculated, calculate!

	this->add(c, lint, lp, l, kappa, tau, kp);
*/
	const Bloch_summed_structure_constant B;
	return B.calc_dot(c, lint, lp, l, kappa, tau, kp, D_m);
}

void Bloch_summed_structure_constant::Container::recalculate_all(const Crystal_t<3, Atom>& c, const lm& lint)
{
    auto old_values = values_m;
    auto old_dot_values = dot_values_m;

    values_m.clear();
    dot_values_m.clear();

	D_m.recalculate_all(c);
    for(const auto& it : old_values){
        auto lp = std::get<0>(it.first);
		auto l = std::get<1>(it.first);
        auto kappa = std::get<2>(it.first);
        auto tau = std::get<3>(it.first);
        auto kp = std::get<4>(it.first);
        this->add(c, lint, lp, l, kappa, tau, kp);
    }

	for(const auto& it : old_dot_values){
        auto lp = std::get<0>(it.first);
		auto l = std::get<1>(it.first);
        auto kappa = std::get<2>(it.first);
        auto tau = std::get<3>(it.first);
        auto kp = std::get<4>(it.first);
        this->dot(c, lint, lp, l, kappa, tau, kp);
    }
}

GSL::Complex Bloch_summed_structure_constant::calc(
	const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa,
	const GSL::Vector& tau,	const GSL::Vector& kp) const
{
	Bloch_sum::Container Ds;
	return this->calc(c, lint, lp, l, kappa, tau, kp, Ds);
}

GSL::Complex Bloch_summed_structure_constant::calc(
	const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa,
	const GSL::Vector& tau,	const GSL::Vector& kp, Bloch_sum::Container& Ds) const
{
	GSL::Complex sum(0., 0.);
	double k_fac = 1;
	if(std::abs(kappa) > 1e-10){
		k_fac = GSL::pow_int(kappa, l.l + lp.l);
		for (int lpp = 0; lpp <= lint.l; lpp++){
			for (int mpp = -lpp; mpp <= lpp; mpp++){
				GSL::Result a = gaunt(l, lp, lm {lpp, mpp});
				if (std::abs(a.val) > 1E-16){
				    sum += GSL::pow_int(-1, lpp)*GSL::pow_int(1./kappa, lpp)*
					a.val*Ds(lm {lpp, mpp}, kappa, c, tau, kp);
				}
			}
		}
	}else{
		int lpp = l.l + lp.l;
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			GSL::Result a = gaunt(l, lp, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-16){
			    sum += GSL::pow_int(-1, lpp)*
				a.val*Ds(lm {lpp, mpp}, kappa, c, tau, kp);
			}
		}
	}

	int sign = l.l % 2 == 0 ? 1 : -1;

	return 4*M_PI*sign*k_fac*sum;
}

GSL::Complex Bloch_summed_structure_constant::calc_dot(
	const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa,
	const GSL::Vector& tau, const GSL::Vector& kp) const
{
	Bloch_sum::Container Ds;
	return this->calc_dot(c, lint, lp, l, kappa, tau, kp, Ds);
}
GSL::Complex Bloch_summed_structure_constant::calc_dot(
	const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa,
	const GSL::Vector& tau, const GSL::Vector& kp, Bloch_sum::Container& Ds) const
{
	GSL::Result a;
	GSL::Complex sum(0., 0.);
	double k_fac = 1;
	if(std::abs(kappa) > 1e-10){
		k_fac = GSL::pow_int(kappa, l.l + lp.l);
		for (int lpp = 0; lpp <= lint.l; lpp++){
			for (int mpp = -lpp; mpp <= lpp; mpp++){
				a = gaunt(l, lp, lm {lpp, mpp});
				if (std::abs(a.val) > 1E-16){
					sum += a.val*
						(2*Ds.dot(lm {lpp,	 mpp}, kappa, c, tau, kp) -
						 (l.l + lp.l - lpp)/GSL::pow_int(kappa, 2)*
						 Ds(lm {lpp, mpp}, kappa, c, tau, kp))*
						GSL::pow_int(-1, lpp)*GSL::pow_int(1./kappa, lpp);
				}
			}
		}
	}else{
		int lpp = l.l + lp.l;
		for (int mpp = -lpp; mpp <= lpp; mpp++){
			a = gaunt(l, lp, lm {lpp, mpp});
			if (std::abs(a.val) > 1E-16){
				sum += a.val*(2*Ds.dot(lm {lpp,	 mpp}, kappa, c, tau, kp));
			}
		}
	}
	int sign = l.l % 2 == 0 ? 1 : -1;

	return 2*M_PI*sign*k_fac*sum;
}
