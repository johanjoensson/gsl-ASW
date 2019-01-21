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
#include "GSLpp/complex_matrix.h"
#include "GSLpp/complex.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/error.h"
#include "augmented_spherical_wave.h"
#include "atomic_quantity.h"
#include "envelope_fun.h"
#include "xc_func.h"
#include "k-mesh.h"

int main()
{
#ifdef DEBUG
std::ofstream numerov_debug;
numerov_debug.open("numerov.debug", std::fstream::out|std::fstream::trunc);
numerov_debug.close();
numerov_debug.open("check_Hankel.dat", std::fstream::out|std::fstream::trunc);
numerov_debug.close();
numerov_debug.open("check_Bessel.dat", std::fstream::out|std::fstream::trunc);
numerov_debug.close();
#endif
	std::cout.precision(12);

	GSL::Error_handler e_handler;
	e_handler.off();


	double kappa = sqrt(0.015);

//	GSL::Vector a = {0.0, 0.5, 0.5}, b = {0.5, 0.0, 0.5}, c = {0.5, 0.5, 0.0};
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

	tau[0] = 0.25;
	tau[1] = 0.25;
	tau[2] = 0.25;
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
	C2.set_Z(6);
	C3.set_Z(6);
	C4.set_Z(6);

	cr.add_atoms(std::vector<Atom> {C1, C2});

	Simulation sim(cr, LDA, kappa);
	sim.set_up_X_matrices();
	K_mesh kmesh(cr.lat.r_lat);
	kmesh.generate_mesh(5, 5, 5);

	Bessel_function j(lm {1, 0});
	std::cout << GSL::pow_int(kappa, -1)*3.3217036494*j(kappa*3.3217036494) << "\n";

	for(GSL::Vector kp : kmesh.k_points){
	// for(GSL::Vector kp : { GSL::Vector {0., 0., 0.}, GSL::Vector { -1.777153, -3.554306, -1.777153 } }){

		std::cout << "k-point " << kp <<
		"\n" << std::string(80, '*') << "\n";
	 	sim.set_up_H(kp);
		sim.set_up_S(kp);
		sim.calc_eigen();
	}

	return 0;
}
