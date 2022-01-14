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
	struct hash<std::tuple<lm, double, double, double, double, double, double, double>>{
		size_t operator()(const std::tuple<lm, double, double, double, double, double, double, double>& key) const
		{
			size_t res = 0;
			res ^= std::hash<lm>()(std::get<0>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<1>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);

			res ^= std::hash<double>()(std::get<2>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<3>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<4>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);

			res ^= std::hash<double>()(std::get<5>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<6>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);
			res ^= std::hash<double>()(std::get<7>(key)) + 0x9e3779b9 + (res<< 6) + (res>> 2);


			return res;
		}
	};
}

class Bloch_sum{

    using Key = std::tuple<lm, double, double, double, double, double, double, double>;
    GSL::Complex calc_d1(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;
    GSL::Complex calc_d2(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;
    GSL::Complex calc_d3(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau) const;
    GSL::Complex calc_d1_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;
    GSL::Complex calc_d2_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;
    GSL::Complex calc_d3_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau) const;
public:

    Bloch_sum() = default;

    GSL::Complex hankel_envelope(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;
    GSL::Complex hankel_envelope_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const;

    GSL::Complex operator()(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const
    {
        return hankel_envelope(l, kappa, c, tau, kp);
    }
    GSL::Complex dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp) const
    {
        return hankel_envelope_dot(l, kappa, c, tau, kp);
    }

	class Container{
	private:
		std::unordered_map<Bloch_sum::Key, GSL::Complex> values_m;
    	std::unordered_map<Bloch_sum::Key, GSL::Complex> dot_values_m;
		void add(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp);
	public:
		Container() /* : values_m(), dot_values_m() */ {}
	    GSL::Complex get(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp);
	    GSL::Complex operator()(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp)
	    {
	        return get(l, kappa, c, tau, kp);
	    }

	    GSL::Complex get_dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp);
	    GSL::Complex dot(const lm l, const double kappa, const Crystal_t<3, Atom>& c, GSL::Vector::Const_View tau, GSL::Vector::Const_View kp)
	    {
	        return get_dot(l, kappa, c, tau, kp);
	    }

		void recalculate_all( const Crystal_t<3, Atom>& c);
	};
};


#endif // BLOCH_SUM_H
