#include <iostream>
#include <fstream>
#include <cmath>
#include "spherical_fun.h"
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
#include "../../GSL-lib/src/complex.h"
#include "../../GSL-lib/src/special_functions.h"

int main()
{
	Logarithmic_mesh mesh(0.5, 500);

	std::cout.precision(12);

	GSL::Vector kp(3), R(3);
    GSL::Complex d2(0., 0.), e(0., 0.);


	GSL::Vector a(3), b(3), c(3);
	a[0] = 0.5; a[1] = 0.5; a[2] = 0.5;
	b[0] = 0.5; b[1] = -0.5; b[2] = 0.5;
	c[0] = 0.5; c[1] = -0.5; c[2] = -0.5;
	Crystal cr(4*a, 4*b, 4*c);
	GSL::Vector tau(3), tau_0(3);
	tau[0] = 0.;
	tau[1] = 0;
	tau[2] = 0.;

	tau_0.copy(tau);
	Atom at(cr.scale*cr.lattice*tau, mesh);
	cr.add_atoms(std::vector<Atom> {at});

	kp[1] = 0.5 ;

	Bloch_sum bs( lm {1, 0}, sqrt(0.015), cr);

	std::cout << "D(";
	std::cout << cr.atoms[0].get_pos() << ", ";
	std::cout << kp << ") = ";
	std::cout << bs.calc_d1(tau, kp) + bs.calc_d2(tau, kp) + bs.calc_d3(tau)
	<< std::endl;

	return 0;
}
