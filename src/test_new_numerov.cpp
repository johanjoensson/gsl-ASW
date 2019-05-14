#include "numerov_solver.h"
#include "schroedinger.h"
#include <vector>
#include <fstream>
#include <GSLpp/special_functions.h>
#include <GSLpp/error.h>
#include <cmath>


int main()
{
	GSL::Error_handler e_handler;
	e_handler.off();
	Numerov_solver sol;
	size_t len = 10000;

	// 1D particle in  box
	int n = 3;
	double l = 1;
	double h_2 = GSL::pow_int(l/(len - 1), 2);
	std::vector<double> res(len, 0), g(len, 0), s(len, 0), left{0., std::sqrt(2/l)*GSL::sin(n*M_PI/l*(-0.5*l + 2./(len - 1))).val},
		right{GSL::pow_int(-1, n - 1)*std::sqrt(2/l)*GSL::sin(n*M_PI/l*(0.5*l - 1./(len -1))).val, 0.};
	for(auto it = g.begin(); it != g.end(); it++){
		*it = GSL::pow_int(n*M_PI/(l), 2)*h_2;
	}
	sol.solve(res.begin(), res.end(), g.begin(), g.end(), s.begin(),
		s.end(), left.begin(), left.end(), right.rbegin(), right.rend());

	std::fstream outfile("numerov_test.dat", std::ios::out);
	for(auto& val : res){
		outfile << val << "\n";
	}
	outfile << "\n\n";

	// 1D harmonic oscillator
	auto fac = [](int n){
		double res = GSL::fact(static_cast<unsigned int>(n)).val;
		return res;};

	l = 6;
	n = 0;
	h_2 = GSL::pow_int(l/(len - 1), 2);
	left = {
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(-l/2, 2)).val*GSL::hermite_phys(n, -l/2).val,
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(-l/2 + 1./(len-1), 2)).val*GSL::hermite_phys(n, -l/2 + 1./(len - 1)).val};
	right = {
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(l/2 - 1./(len - 1), 2)).val* GSL::hermite_phys(n, l/2 - 1./(len - l)).val,
		1./std::sqrt(GSL::pow_int(2, n)*fac(n))*1./std::sqrt(std::sqrt(M_PI))*
			GSL::exp(-0.5*GSL::pow_int(l/2, 2)).val*GSL::hermite_phys(n, l/2).val
	 };

	int i = 0;
	for(auto it = g.begin(); it != g.end(); it++, i++){
		*it = (2*n + 1 - GSL::pow_int(
			i*l/(len - 1) - l/2, 2))*h_2;
	}
	sol.solve(res.begin(), res.end(), g.begin(), g.end(), s.begin(),
		s.end(), left.begin(), left.end(), right.rbegin(), right.rend());

	for(auto& val : res){
		outfile << val << "\n";
	}
	outfile << "\n\n";

	return 0;
}
