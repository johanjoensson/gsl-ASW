#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include "gaunt.h"
#include "numerov_solver.h"

int l;
int z;
double kappa = std::sqrt(0.015);

double v_at(double r)
{
	return -2.*z/r;
}

double v_eff(double r)
{
	return l*(l+1.)/(r*r) + v_at(r);
}



int main()
{
	z = 8;
	l = 2;
	int n = 4;

	unsigned int len = 1000;
	Logarithmic_mesh mesh(10.0, len);
	Numerov_solver sol;
	sol.set_v_eff(&v_eff);
	sol.set_v_at(&v_at);

	// Initial guess of energy eigenvalue (theoretical value = -z^2/n^2)
	double energy = -z*z*1./(n*n) -2.;

	std::vector<double> res(len,0);

	std::vector<double> icond_r { 
		  gsl_pow_int(kappa, -l)*gsl_sf_bessel_jl(l, kappa*mesh.r[len-1])/std::sqrt(mesh.drx[len-1]),
		  gsl_pow_int(kappa, -l)*gsl_sf_bessel_jl(l, kappa*mesh.r[len-2])/std::sqrt(mesh.drx[len-2])
	}; 

	std::vector<double> icond_l = { gsl_pow_int(mesh.r[0], l+1)/std::sqrt(mesh.drx[0]), 
		  gsl_pow_int(mesh.r[1], l+1)/std::sqrt(mesh.drx[1]), 
		  gsl_pow_int(mesh.r[2], l+1)/std::sqrt(mesh.drx[2]) };

	res = sol.solve(mesh, icond_l, icond_r, energy, n - l - 1);
	std::cout << "Returned energy = " << energy << std::endl;

	std::ofstream file;
	file.open("WaveFun.dat");
	for(unsigned int i = 0; i < res.size(); i++ ){
		file << mesh.r[i] << " " << res[i] << std::endl;
	}
	file.close();

	return 0;
}
