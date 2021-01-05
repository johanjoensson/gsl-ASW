#include <iostream>
#include <algorithm>
#include <tuple>
#include <cmath>
#include <gaunt.h>
#include <spherical_fun.h>
#include <GSLpp/special_functions.h>
#include <GSLpp/complex.h>

#include <chrono>

/*******************************************************************************
* Check triangle identity for l-quantum numbers l1.l, l2.l and l3.l            *
* | l1.l - l2.l | < l3.l < | l1.l + l2.l |
*******************************************************************************/
constexpr inline bool triangle(lm l1, lm l2, lm l3)
{
	return (l3.l >= std::abs(l1.l - l2.l) && l3.l <= std::abs(l1.l + l2.l));
}

/*******************************************************************************
* Get coefficients and lm pair for complex spherical harmonics that make up    *
* the real cubic harmonic with lm pair l.                                      *
*******************************************************************************/
std::vector<std::pair<GSL::Complex, lm>> get_complex_sph (const lm& l)
{
	std::vector<std::pair<GSL::Complex, lm>> res;
	int l1 = l.l, l2 = l.l;
	int m1 = l.m, m2 = l.m;
	GSL::Complex c1 = 0, c2 = 0;

	if(m1 < 0){
		m2 = -m2;
		c1 = GSL::Complex(0, 1)/sqrt(2);
		c2 = GSL::pow_int(-1, m2 + 1)*c1;
		res.push_back({c1, {l1, m1}});
		res.push_back({c2, {l2, m2}});
	}else if (m1 > 0){
		m1 = -m1;
		c1 = 1/sqrt(2);
		c2 = GSL::pow_int(-1, m2)*c1;
		res.push_back({c1, {l1, m1}});
		res.push_back({c2, {l2, m2}});
	}else{
		c1 = 1;
		c2 = 0;
		res.push_back({c1, {l1, m1}});
	}
	return res;
}


inline GSL::Result gaunt(lm l1, lm l2, lm l3)
{
	double C = std::sqrt((2.*l1.l + 1)*(2.*l2.l + 1)*(2.*l3.l + 1)/(4*M_PI));
	return  C*GSL::wigner_3j(l1.l, l2.l, l3.l,  0,    0,    0  )*
	  			GSL::wigner_3j(l1.l, l2.l, l3.l, l1.m, l2.m, l3.m);
}

/*******************************************************************************
* NAIVE IMPLEMENTATION, BUT IT SEEMS TO WORK                                   *
*******************************************************************************/
GSL::Result real_gaunt(lm l1, lm l2, lm l3)
{
	GSL::Complex res(0, 0);
	double err = 0;
	for(auto c1 : get_complex_sph(l1)){
		for(auto c2 : get_complex_sph(l2)){
			for(auto c3 : get_complex_sph(l3)){
				GSL::Result g = gaunt(c1.second, c2.second, c3.second);
				res += c1.first*c2.first*c3.first*g.val;
				err += g.err;
			}
		}
	}
	if(abs(res.im()) > 1e-6){
		std::cerr << "Real Gaunt coefficient is " << res << ", with a rather large imaginary part!\n";
	}
	return GSL::Result(res.re(), err);
}
