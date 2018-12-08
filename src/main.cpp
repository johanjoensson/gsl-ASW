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

	GSL::Vector a = {0.0, 0.5, 0.5}, b = {0.5, 0.0, 0.5}, c = {0.5, 0.5, 0.0};
	std::cout << a << std::endl;
	std::cout << b << std::endl;
	std::cout << c << std::endl;
	Crystal cr(6*a, 6*b, 6*c);

	std::cout << cr.lat.scale*cr.lat.lat << std::endl;

	std::cout << "Calculating Rmax" << std::endl;
	double Rmax = cr.calc_Rmax(1e-14, sqrt(0.015), lm {4, 0});
	std::cout << "Calculating Kmax" << std::endl;
	double Kmax = cr.calc_Kmax(1e-14, sqrt(0.015), lm {4, 0});


	std::cout << "Setting R-vectors" << std::endl;
	cr.set_Rn(Rmax);
	std::cout << "Setting K-vectors" << std::endl;
	cr.set_Kn(Kmax);

	std::cout << "Setting up atoms" << std::endl;
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

//	cr.add_atoms(std::vector<Atom> {C1, C2, C3, C4});
	cr.add_atoms(std::vector<Atom> {C1, C2});

	Simulation sim(cr, LDA, sqrt(0.015));
	sim.set_up_X_matrices();
	K_mesh kmesh(cr.lat.r_lat);
	kmesh.generate_mesh(2, 2, 2);
	for(GSL::Vector kp : {GSL::Vector(3)}){
		std::cout << kp << std::endl;
		sim.set_up_H(kp);
		sim.set_up_S(kp);
		sim.calc_eigen();
	}

	return 0;
}
