#ifndef BLOCH_SUM_H
#define BLOCH_SUM_H
#include "spherical_fun.h"
#include "crystal.h"
#include "atom.h"
#include "GSLpp/vector.h"
#include "GSLpp/complex.h"
#include <unordered_map>
#include <tuple>
#include <functional>


namespace std{
	template<>
	struct hash<std::tuple<lm, double, GSL::Vector, GSL::Vector>>{
		size_t operator()(const std::tuple<lm, double, GSL::Vector, GSL::Vector>& key) const
		{
			size_t res = 0;
			res ^= std::hash<lm>()(std::get<0>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<1>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<GSL::Vector>()(std::get<2>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<GSL::Vector>()(std::get<3>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			return res;
		}
	};
}

class Bloch_sum{

    using Key = std::tuple<lm, double, GSL::Vector, GSL::Vector>;
    static std::unordered_map<Key, GSL::Complex> values_m;
    static std::unordered_map<Key, GSL::Complex> dot_values_m;

    GSL::Complex calc_d1(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;
    GSL::Complex calc_d2(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;
    GSL::Complex calc_d3(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau) const;
    GSL::Complex calc_d1_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;
    GSL::Complex calc_d2_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;
    GSL::Complex calc_d3_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau) const;
public:
    Bloch_sum(){}

    GSL::Complex hankel_envelope(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;
    GSL::Complex hankel_envelope_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const;

    GSL::Complex operator()(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
    {
        return hankel_envelope(l, kappa, c, tau, kp);
    }
    GSL::Complex dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, const GSL::Vector& tau, const GSL::Vector& kp) const
    {
        return hankel_envelope_dot(l, kappa, c, tau, kp);
    }
};


#endif // BLOCH_SUM_H
