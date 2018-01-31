#include <iostream>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "numerov_solver.h"

int l = 2;

double real_spherical_hankel(int l, double x)
{
  double exp = gsl_sf_exp(-x);
  double k = gsl_sf_bessel_kl_scaled(l, x);

  return exp*k;
}

double empty_pot(double r)
{
	return 0.*r;
}

/* To change later */
double v_eff(double r)
{
	return l*(l+1.)/(r*r) - 2./r;
}

Numerov_solver::Numerov_solver()
	: v()
{
	this->v_ext = &empty_pot;
}

void Numerov_solver::set_v_ext(double (*new_v_ext)(double r))
{
	this->v_ext = new_v_ext;
}

std::vector<double> Numerov_solver::solve_left(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E)
{
	uint l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;

	for(uint i = res.size() - 1; i > res.size() - 1 - l; i--){
		res[i] = init_cond[res.size() - 1 - i];
	}

	for(int i = res.size() - 1 - l; i >= i_lr; i--){
		g = mesh.drx[i]*mesh.drx[i]*(E - v_eff(mesh.r[i])) - mesh.A*mesh.A/4.;
		g1 = mesh.drx[i+1]*mesh.drx[i+1]*(E - v_eff(mesh.r[i+1])) - mesh.A*mesh.A/4.;
		g2 = mesh.drx[i+2]*mesh.drx[i+2]*(E - v_eff(mesh.r[i+2])) - mesh.A*mesh.A/4.;
		res[i] = (2*res[i+1]*(1 - 5*g1/12) - res[i+2]*(1 + g2/12))/(1 + g/12);
	}

	return res;
}

std::vector<double> Numerov_solver::solve_right(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E)
{
	uint l = init_cond.size();
	std::vector<double> res(mesh.r.size(),0);
	double g, g1, g2;

	for(uint i = 0; i < l; i++){
		res[i] = init_cond[i];
	}

	for(int i = l; i <= i_lr; i++){
		g = mesh.drx[i]*mesh.drx[i]*(E - v_eff(mesh.r[i])) - mesh.A*mesh.A/4.;
		g1 = mesh.drx[i-1]*mesh.drx[i-1]*(E - v_eff(mesh.r[i-1])) - mesh.A*mesh.A/4.;
		g2 = mesh.drx[i-2]*mesh.drx[i-2]*(E - v_eff(mesh.r[i-2])) - mesh.A*mesh.A/4.;
		res[i] = (2*res[i-1]*(1 - 5*g1/12) - res[i-2]*(1 + g2/12))/(1 + g/12);
	}

	return res;
}

int main()
{
	uint len = 1000;
	Logarithmic_mesh mesh(20.0, len);
	Numerov_solver sol;

	std::vector<double> res(len,0);

	std::vector<double> icond { 
		  gsl_pow_int(sqrt(0.015), l+1)*real_spherical_hankel(l, sqrt(0.015)*mesh.r[len-1])/sqrt(mesh.drx[len-1]),
		  gsl_pow_int(sqrt(0.015), l+1)*real_spherical_hankel(l, sqrt(0.015)*mesh.r[len-2])/sqrt(mesh.drx[len-2]) }; 

	std::vector<double> left = sol.solve_left(mesh, len-25, icond, -0.112);

	icond = { gsl_pow_int(mesh.r[0], l+1)/sqrt(mesh.drx[0]), 
		  gsl_pow_int(mesh.r[1], l+1)/sqrt(mesh.drx[1]), 
		  gsl_pow_int(mesh.r[2], l+1)/sqrt(mesh.drx[2]) };

	std::vector<double> right = sol.solve_right(mesh, len-25, icond, -0.112);

	double scale = right[len-25]/left[len-25];
	left[len-25] *= 0.5; /* Avoid double counting */
	right[len-25] *= 0.5; /* Avoid double counting */
	for(uint i = 0; i < len; i++ ){
		res[i] = scale*left[i] + right[i];
	}

	for(uint i = 0; i < len; i++)
 		std::cout << mesh.r[i] << " " << res[i]*sqrt(mesh.drx[i])/mesh.r[i] << std::endl;


	return 0;
}
