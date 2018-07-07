#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
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
//	gsl_set_error_handler_off();
	double vol = 90;
	double eta = exp_gsl(gsl_sf_log(6.5) + 2./3*gsl_sf_log(4*M_PI/3) - 2./3*gsl_sf_log(vol));
	Logarithmic_mesh mesh(0.5, 500);
	std::vector<double> ew(mesh.r.size(), 0.);
	std::vector<double> ew_comp(mesh.r.size(), 0.);
	Ewald_integral I;
	I.set_kappa(sqrt(0.015));
	I.set_ewald_param(eta);

	gsl_vector *tau, *kp, *R, *a1, *a2, *a3, *tau_orig;
    gsl_complex d2 = gsl_complex_rect(0., 0.), e = gsl_complex_rect(0., 0.);

	tau_orig = gsl_vector_alloc(3);
	tau = gsl_vector_alloc(3);
	kp = gsl_vector_alloc(3);
	R = gsl_vector_alloc(3);
	a1 = gsl_vector_alloc(3);
	a2 = gsl_vector_alloc(3);
	a3 = gsl_vector_alloc(3);


	double x = 0., y = 0., z = 0, r = 0.;

    double dot = 0.;
	double tmp = 0.;
	int l = 5, m = -3;

	gsl_vector_set(kp, 0, 0.5);
	gsl_vector_set(kp, 1, 0.5);
	gsl_vector_set(kp, 2, 0.5);

	gsl_vector_set(tau, 0, 0.2);
	gsl_vector_set(tau, 1, 0.86);
	gsl_vector_set(tau, 2, 0.1);
	gsl_vector_memcpy(tau_orig, tau);
	for(int n = 1; n < 100; n++){
	d2 = gsl_complex_rect(0., 0.);
	for(int i = 1; i < n; i++){
		for(int j = 1; j < n; j++){
			for(int k = 1; k < n; k++){
				gsl_vector_memcpy(tau, tau_orig);
				gsl_vector_set_zero(R);
				gsl_vector_set_basis(a1, 0);
				gsl_vector_set_basis(a2, 1);
				gsl_vector_set_basis(a3, 2);

				gsl_vector_scale(a1, double(4*i));
				gsl_vector_scale(a2, double(4*j));
				gsl_vector_scale(a3, double(5*k));

/*				R = *a1 + *a2;
				gsl_vector_ass(R, a3);
*/
				gsl_vector_add(R, a1);
				gsl_vector_add(R, a2);
				gsl_vector_add(R, a3);
				gsl_vector_sub(tau, R);
				x = gsl_vector_get(tau, 0);
				y = gsl_vector_get(tau, 1);
				z = gsl_vector_get(tau, 2);
				r = sqrt(x*x + y*y + z*z);


				gsl_blas_ddot(R, kp, &dot);
				tmp = gsl_pow_int(r, l)*cubic_harmonic(lm {l, m}, *tau)*I.ewald_int(lm {l, m}, r);
				e = gsl_complex_exp(gsl_complex_rect(0., dot));
				d2 = gsl_complex_add(d2, gsl_complex_mul_real(e, tmp));
			}
		}
	}
	d2 = gsl_complex_mul_real(d2, 2./M_SQRTPI*gsl_pow_int(2, l));

	std::cout << "N = " << n << ": D2([";
	std::cout << gsl_vector_get(tau_orig, 0) << ", " << gsl_vector_get(tau_orig, 1) << ", " << gsl_vector_get(tau_orig, 2) << "], [";
	std::cout << gsl_vector_get(kp, 0) << ", " << gsl_vector_get(kp, 1) << ", " << gsl_vector_get(kp, 2) << "]) = ";
	std::cout << d2<< std::endl;
    }
/*	ew = I.evaluate(lm {l, 0}, mesh);
	ew_comp = I.evaluate_comp(lm {l, 0}, mesh);
	double r = 0, ir = 0, ir2 = 0, ir3 = 0, kappa = sqrt(0.015), kappa_fac = gsl_pow_int(kappa, l+1);
	double two_l = gsl_pow_int(2., l), r_l = 0;
	for(unsigned int i = 1; i < mesh.r.size(); i++){
		r = mesh.r[i];
		ir = 1./r;
		ir2 = gsl_pow_int(ir, 2);
		ir3 = gsl_pow_int(ir, 3);
		r_l = gsl_pow_int(r, l);

		std::cout << r << " " << two_l*r_l*2./M_SQRTPI*(ew[i] + ew_comp[i]) << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) << std::endl;
		std::cout << r << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) - exp_gsl(-kappa*r)/r << std::endl;
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
