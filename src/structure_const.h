#ifndef STRUCTURE_CONST_H
#define STRUCTURE_CONST_H
#include <iostream>
#include <functional>
#include <unordered_map>
#include <tuple>
//#include <gsl/gsl_vector.h>
#include "utils.h"
#include "bloch_sum.h"
#include "GSLpp/vector.h"
#include "GSLpp/complex.h"

/***************************************************************************//**
* A class for representing structure constants\n
* Contains:\n
* __l_int__ - Maximum orbital angular momentum to be included for Bessel
* expansions\n __l_low__ - Maximumm orbital angular momentum to be included for
* Hankel expansions\n __l1__, __l2__ - Orbital angular momenta to couple via the
* structure constant\n
* __kappa__ - Energy parameter used (usually kappa^2 = -0.015)\n
* __r__ - Position of atom\n
* __val__ - Value of the structure constant\n
* __dk_val__ - Value of energy derivative of the structure constant\n
*******************************************************************************/
class Structure_constant{
		int l_int;
	    lm l1, l2;
		double kappa;
	public:
		Structure_constant(int l_int_n, double kappa_n, lm l1_n, lm l2_n);
		Structure_constant(int l_int_n, lm l1_n, lm l2_n);
		Structure_constant(lm l1_n, lm l2_n);
		// Structure_constant() = default;

		double operator()(const GSL::Vector& r) const;
		double dot(const GSL::Vector& r) const;
		friend std::ostream& operator << ( std::ostream&,
			const Structure_constant& );
};
namespace std{
	template<>
	struct hash<std::tuple<lm, lm, double, GSL::Vector, GSL::Vector>>{
		size_t operator()(const std::tuple<lm, lm, double, GSL::Vector, GSL::Vector>& key) const
		{
			size_t res = 0;
			res ^= std::hash<lm>()(std::get<0>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<lm>()(std::get<1>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<2>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<GSL::Vector>()(std::get<3>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<GSL::Vector>()(std::get<4>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			return res;
		}
	};
}

class Bloch_summed_structure_constant{
	using Key = std::tuple<lm, lm, double, GSL::Vector, GSL::Vector>;
public:
	Bloch_summed_structure_constant(){}
	GSL::Complex calc(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp, Bloch_sum::Container& Ds) const;
	GSL::Complex calc(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp) const;
	GSL::Complex operator()(const Crystal_t<3, Atom>& c, const lm lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp, Bloch_sum::Container& Ds) const
	{
		return this->calc(c, lint, lp, l, kappa, tau, kp, Ds);
	}
	GSL::Complex operator()(const Crystal_t<3, Atom>& c, const lm lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp) const
	{
		return this->calc(c, lint, lp, l, kappa, tau, kp);
	}

	GSL::Complex calc_dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp, Bloch_sum::Container& Ds) const;
	GSL::Complex calc_dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp) const;
	GSL::Complex dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp, Bloch_sum::Container& Ds) const
	{
		return this->calc_dot(c, lint, lp, l, kappa, tau, kp, Ds);
	}
	GSL::Complex dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp) const
	{
		return this->calc_dot(c, lint, lp, l, kappa, tau, kp);
	}

	class Container {
	private:
		std::unordered_map<Bloch_summed_structure_constant::Key, GSL::Complex> values_m;
		std::unordered_map<Bloch_summed_structure_constant::Key, GSL::Complex> dot_values_m;
		Bloch_sum::Container D_m;

		void add(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp);
	public:
		Container() : values_m(), dot_values_m(), D_m() {}
		GSL::Complex get(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp);
		GSL::Complex get_dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp);
		GSL::Complex operator()(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp)
		{
			return this->get(c, lint, lp, l, kappa, tau, kp);
		}
		GSL::Complex dot(const Crystal_t<3, Atom>& c, const lm& lint, const lm& lp, const lm& l, const double kappa, const GSL::Vector& tau, const GSL::Vector& kp)
		{
			return this->get_dot(c, lint, lp, l, kappa, tau, kp);
		}

		void recalculate_all(const Crystal_t<3, Atom>& c, const lm& lint);
	};
};

std::ostream& operator << ( std::ostream&, const Structure_constant& );
#endif //STRUCTURE_CONST_H
