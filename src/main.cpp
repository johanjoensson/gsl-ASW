#include <iostream>
#include <fstream>
#include <cmath>
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include "numerov_solver.h"
#include "atom.h"
#include "crystal.h"
#include "ewald_int.h"
#include "gaunt.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/complex.h"
#include "../../GSL-lib/src/special_functions.h"


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
	double vol = 90;
	double eta = GSL::exp(GSL::log(6.5) + 2./3*GSL::log(4*M_PI/3) -
	2./3*GSL::log(vol)).val;
	Logarithmic_mesh mesh(0.5, 500);
	std::vector<double> ew(mesh.r.size(), 0.);
	std::vector<double> ew_comp(mesh.r.size(), 0.);
	Ewald_integral I;
	I.set_kappa(sqrt(0.015));
	I.set_ewald_param(eta);

	std::cout.precision(12);

	GSL::Vector tau, kp, R, a1, a2, a3, tau_orig;
    GSL::Complex d2(0., 0.), e(0., 0.);

	tau_orig = GSL::Vector(3);
	tau = GSL::Vector(3);
	kp = GSL::Vector(3);
	R = GSL::Vector(3);
	a1 = GSL::Vector(3);
	a2 = GSL::Vector(3);
	a3 = GSL::Vector(3);

	a1[0] = 0.5;
	a1[1] = 0.5;
	a1[2] = 0;

	a2[0] = 0;
	a2[1] = 0.5;
	a2[2] = 0.5;

	a3[0] = 0.5;
	a3[1] = 0.5;
	a3[2] = 0;

    double dot = 0.;
	double r = 0;
	GSL::Result tmp;
	int l = 3, m = -2;

	kp[0] = 0.;
	kp[1] = 0.;
	kp[2] = 0.;

	tau[0] = 0.;
	tau[1] = 0.;
	tau[2] = 0.;


	tau_orig.copy(tau);
	for(int n = 1; n < 20; n++){
	d2 = GSL::Complex(0., 0.);
	for(int i = 1; i < n; i++){
		for(int j = 1; j < n; j++){
			for(int k = 1; k < n; k++){
				R  = i*a1 + j*a2 + k*a3;
				tau = tau_orig + R;
				r = tau.norm();


				dot = GSL::dot(tau, R);
				tmp = GSL::pow_int(r, l)*cubic_harmonic(lm {l, m}, tau)*I.ewald_int(lm {l, m}, r);
				e = GSL::exp(GSL::Complex(0., dot));
				d2 += e*tmp.val;
			}
		}
	}
	d2 *= 2./M_SQRTPI*GSL::pow_int(2, l);

	std::cout << "N = " << n << ": D2(";
	std::cout << tau_orig << ", ";
	std::cout << kp << ") = ";
	std::cout << d2 << std::endl;
    }
	std::cout << std::endl;

	ew = I.evaluate(lm {l, 0}, mesh);
	ew_comp = I.evaluate_comp(lm {l, 0}, mesh);
//	double kappa = sqrt(0.015), kappa_fac = GSL::pow_int(kappa, l+1);
//	double two_l = GSL::pow_int(2., l), r_l = 0;
	for(unsigned int i = 1; i < mesh.r.size(); i++){
//		r = mesh.r[i];
//		r_l = GSL::pow_int(r, l);

//		std::cout << r << " " << two_l*r_l*2./M_SQRTPI*(ew[i] + ew_comp[i]) << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) << std::endl;
//		std::cout << r << " " << kappa_fac*real_spherical_hankel(lm {l, 0}, kappa*r) - GSL::exp(-kappa*r)/r << std::endl;
	}
	/*
	GSL::Vector a(3);
	double x = 0., y = 0., z = 0.;
	for (int theta = 0; theta <= 180; theta += 10){
		for(int phi = -180; phi <= 180; phi += 10){
			x = GSL::sin(theta*M_PI/180).val*GSL::cos(phi*M_PI/180).val;
			y = GSL::sin(theta*M_PI/180).val*GSL::sin(phi*M_PI/180).val;
			z = GSL::cos(theta*M_PI/180).val;
			a[0] = x;
			a[1] = y;
			a[2] = z;
			std::cout << theta*M_PI/180 <<" " << phi*M_PI/180 << " " << cubic_harmonic(lm {2, 0}, a) << std::endl;
		}
	}
	*/
	return 0;
}
