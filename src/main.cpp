#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
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
	a[0] = -0.5; a[1] = 0.5; a[2] = 0.5;
	b[0] = 0.5; b[1] = -0.5; b[2] = 0.5;
	c[0] = 0.5; c[1] = 0.5; c[2] = -0.5;
	Crystal cr(1*a, 1*b, 1*c);
	size_t Nk = cr.calc_nk(1e-6, sqrt(0.015), lm {3, 0});
	size_t Nr = cr.calc_nr(1e-6, sqrt(0.015), lm {3, 0});
	std::cout << "Nr = " << Nr << std::endl;
	cr.calc_Rn(Nr/3);
	std::cout << "Nk = " << Nk << std::endl;
	cr.calc_Kn(Nk/3);


	GSL::Vector tau(3);
	tau[0] = 0.5;
	tau[1] = 0.5;
	tau[2] = 0.5;
	kp[0] = 0.;

	Atom at(cr.lat.scale*cr.lat.lat*tau);
	cr.add_atoms(std::vector<Atom> {at});

	Bloch_sum bs( lm {3, 2}, sqrt(0.015), cr);

	Bloch_summed_structure_constant B(cr, lm {3, 2}, lm {1, -1});

	std::vector<std::vector<Atom>> nn = cr.calc_nearest_neighbours();
	cr.atoms[0].set_MT(nn[0][0].get_pos().norm()/2);
	for(size_t i = 0; i < nn.size(); i++){
		std::cout << "Atom " << i << " : " << std::endl;
		for(size_t j = 0; j < nn[i].size(); j++){
			std::cout << "-> " << j + i << " dist = " << nn[i][j].get_pos().norm() << ";" << std::endl;
		}
		std::cout << std::endl;
	}

	return 0;
}
