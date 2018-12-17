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
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/matrix.h"
#include "../../GSL-lib/src/complex_matrix.h"
#include "../../GSL-lib/src/complex.h"
#include "../../GSL-lib/src/special_functions.h"
#include "../../GSL-lib/src/error.h"
#include "augmented_spherical_wave.h"
#include "atomic_quantity.h"
#include "envelope_fun.h"
#include "xc_func.h"
#include "k-mesh.h"

int main()
{
#ifdef DEBUG
std::ofstream numerov_debug;
numerov_debug.open("numerov.debug", std::fstream::out);
numerov_debug.close();
#endif
	std::cout.precision(12);

	GSL::Error_handler e_handler;
	e_handler.off();

	//GSL::Vector a = {0.0, 0.5, 0.5}, b = {0.5, 0.0, 0.5}, c = {0.5, 0.5, 0.0};
	GSL::Vector a = {1.0, 0.0, 0.0}, b = {0.0, 1.0, 0.0}, c = {0.0, 0.0, 1.0};
	std::cout << a << std::endl;
	std::cout << b << std::endl;
	std::cout << c << std::endl;
	Crystal cr(14*a, 14*b, 14*c);

	std::cout << cr.lat.scale*cr.lat.lat << std::endl;

	std::cout << "Calculating Rmax" << std::endl;
	double Rmax = cr.calc_Rmax(1e-12, sqrt(0.015), lm {4, 0});
	std::cout << "Calculating Kmax" << std::endl;
	double Kmax = cr.calc_Kmax(1e-12, sqrt(0.015), lm {4, 0});


	std::cout << "Setting R-vectors" << std::endl;
	cr.set_Rn(Rmax);
	std::cout << "Setting K-vectors" << std::endl;
	cr.set_Kn(Kmax);

	std::cout << "Setting up atoms" << std::endl;
	GSL::Vector tau(3);
	Atom C1;
	C1.set_pos(tau*cr.lat.scale*cr.lat.lat);

	tau[0] = 0.5;
	tau[1] = 0.5;
	tau[2] = 0.;
	Atom C2;
	C2.set_pos(tau*cr.lat.scale*cr.lat.lat);

	tau[0] = 0.;
	tau[1] = 0.5;
	tau[2] = 0.5;
	Atom C3;
	C3.set_pos(tau*cr.lat.scale*cr.lat.lat);

	tau[0] = 0.5;
	tau[1] = 0.;
	tau[2] = 0.5;
	Atom C4;
	C4.set_pos(tau*cr.lat.scale*cr.lat.lat);

	C1.set_Z(6);
	C2.set_Z(1);
	C3.set_Z(6);
	C4.set_Z(6);

	cr.add_atoms(std::vector<Atom> {C1});

	Simulation sim(cr, LDA, sqrt(0.015));
	sim.set_up_X_matrices();
	K_mesh kmesh(cr.lat.r_lat);
	kmesh.generate_mesh(4, 4, 4);

	// for(GSL::Vector kp : kmesh.k_points){
	for(GSL::Vector kp : { GSL::Vector {0., 0., 0.}, GSL::Vector {-3.926991, -0.785398, -0.785398}}){

		std::cout << "k-point " << kp << std::endl;
		// sim.set_up_H(kp);
		// sim.set_up_S(kp);
//		sim.calc_eigen();
	}

	Ewald_integral I;
	I.set_kappa(sqrt(0.015));
	I.set_ewald_param(16.5);

	Envelope_Hankel h(C1, lm {0, 0}, sqrt(0.015));

	std::cout << I.ewald_int(lm {0, 0}, 0.5) + 2./M_SQRTPI*I.comp_ewald_int(lm {0, 0}, 0.5) << " " << h.barred_fun(0.5) << std::endl;
	return 0;
}
