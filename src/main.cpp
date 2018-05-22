#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include "numerov_solver.h"
#include "atom.h"
#include "crystal.h"
#include "ewald_int.h"

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
	double vol = 1;
	double eta = gsl_sf_exp(gsl_sf_log(6.5) + 2./3*gsl_sf_log(4*M_PI/3) - 2./3*gsl_sf_log(vol));
	Logarithmic_mesh mesh(0.5, 500);
	std::vector<double> ew(mesh.r.size(), 0.);
	std::vector<double> ew_comp(mesh.r.size(), 0.);
	Ewald_integral I;
	I.set_kappa(sqrt(0.015));
	I.set_ewald_param(eta);

	gsl_vector R, a1, a2, a3;

	R = gsl_vector_alloc(3);
	a1 = gsl_vector_alloc(3);
	a3 = gsl_vector_alloc(3);
	a2 = gsl_vector_alloc(3);

	gsl_vector_set_basis(a1, 0);
	gsl_vector_set_basis(a2, 1);
	gsl_vector_set_basis(a3, 2);
	gsl_vector_set_zero(R);

    double d2 = 0.;

	int l = 2;

	for(int i = 0; i < 10; i++){
		for(int j = 0; j < 10; j++){
			for(int k = 0; k< 10; k++){
				gsl_vector_add(R, gsl_vector_scale(a1, double(i)));
				gsl_vector_add(R, gsl_vector_scale(a2, double(j)));
				gsl_vector_add(R, gsl_vector_scale(a3, double(k)));
			}
		}
	}
	ew = I.evaluate(lm {l, 0}, mesh);
	ew_comp = I.evaluate_comp(lm {l, 0}, mesh);
/*	double r = 0, ir = 0, ir2 = 0, ir3 = 0, kappa = sqrt(0.015), kappa_fac = gsl_pow_int(kappa, l+1);
	double two_l = gsl_pow_int(2., l), r_l = 0;
	for(unsigned int i = 1; i < mesh.r.size(); i++){
		r = mesh.r[i];
		ir = 1./r;
		ir2 = gsl_pow_int(ir, 2);
		ir3 = gsl_pow_int(ir, 3);
		r_l = gsl_pow_int(r, l);

		std::cout << r << " " << two_l*r_l*2./M_SQRTPI*(ew[i] + ew_comp[i]) << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) << std::endl;
		std::cout << r << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) - gsl_sf_exp(-kappa*r)/r << std::endl;
	}
	gsl_vector *a = gsl_vector_alloc(3);
	double x = 0, y = 0, z = 0;
	for (double theta = 0.; theta <= M_PI; theta += 0.05){
		for(double phi = -M_PI; phi <= M_PI; phi += 0.05){
				x = sin(theta)*cos(phi);
				y = sin(theta)*sin(phi);
				z = cos(theta);
				gsl_vector_set(a, 0, x);
				gsl_vector_set(a, 1, y);
				gsl_vector_set(a, 2, z);
				std::cout << theta <<" " << phi << " " << cubic_harmonic(lm {2, 0}, *a) << std::endl;
		}
	}
	*/
	return 0;
}
