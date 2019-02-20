#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include "simulation.h"
#include "spherical_fun.h"
#include "utils.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include "numerov_solver.h"
#include "atom.h"
#include "crystal.h"
#include "bloch_sum.h"
#include "ewald_int.h"
#include "gaunt.h"
#include "GSLpp/vector.h"
#include "GSLpp/matrix.h"
#include "GSLpp/complex.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/error.h"
#include "augmented_spherical_wave.h"
#include "atomic_quantity.h"
#include "envelope_fun.h"
#include "xc_func.h"
#include "k-mesh.h"

#include <gsl/gsl_sf_erf.h>

int main()
{
	std::cout.precision(12);

	GSL::Error_handler e_handler;
	e_handler.off();


	Hankel_function h;
	Integral_Hankel_function hi;

	for(int l = 0; l < 3; l++){
		for(int m = -l; m <= l; m++){
			h.set_l(lm{l,m});
			hi.set_l(lm{l,m});
			for(double x = 0.001; x < 9; x += 0.1){
				std::cout << x << " " << h(x) << " " << hi(x) << " " << 1 - hi(x)/h(x) <<"\n";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
	std::cout << "\n";
	double kappa = std::sqrt(0.015);

	GSL::Vector a = {1.0, 0.0, 0.0}, b = {0.0, 1.0, 0.0}, c = {0.0, 0.0, 1.0};
	std::cout << a << "\n";
	std::cout << b << "\n";
	std::cout << c << "\n";
	Crystal cr(4*a, 4*b, 4*c);

	std::cout << cr.lat.scale*cr.lat.lat << "\n";

	std::cout << "Calculating Rmax" << "\n";
	double Rmax = cr.calc_Rmax(5e-14, kappa, lm {5, 0});
	std::cout << "Calculating Kmax" << "\n";
	double Kmax = cr.calc_Kmax(5e-14, kappa, lm {5, 0});
	std::cout << "Rmax = " << Rmax << ", Kmax = " << Kmax << "\n";


	std::cout << "Setting R-vectors" << "\n";
	cr.set_Rn(Rmax);
	std::cout << "Setting K-vectors" << "\n";
	cr.set_Kn(Kmax);

	std::cout << "Setting up atoms" << "\n";
	GSL::Vector tau(3);
	Atom C1;
	C1.set_pos(tau*cr.lat.scale*cr.lat.lat);
	C1.set_Z(6);
	cr.add_atoms(std::vector<Atom> {C1});

	Simulation sim(cr, NONE, kappa);
	sim.set_up_X_matrices();
	K_mesh kmesh(cr.lat.r_lat);

	std::vector<GSL::Vector> k_path = {
		{0.,0.,0.},			// Gamma
		{0.5, 0.5, 0.5},	// R
		{0.5, 0.5, 0.},		// M
		{0., 0.5, 0.}, 		// X
		{0., 0., 0.}		// Gamma
	};

	kmesh.generate_mesh(k_path, 50);

	std::pair<GSL::Matrix_cx, GSL::Vector> eigen;
	for(const GSL::Vector& kp : kmesh.k_points){
		std::cout << "k-point " << kp <<
		"\n" << std::string(80, '*') << "\n";
	 	sim.set_up_H(kp);
		sim.set_up_S(kp);
		eigen = sim.calc_eigen();
		sim.add_eigvals(kp, eigen.second);
	}
	sim.print_eigvals(kmesh);

	return 0;
}
