#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
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
#include "../../GSL-lib/src/complex.h"
#include "../../GSL-lib/src/special_functions.h"
#include "augmented_spherical_wave.h"
#include "atomic_quantity.h"

int main()
{
	Logarithmic_mesh mesh(0.5, 500);

	std::cout.precision(12);

	GSL::Vector kp(3), R(3);
    GSL::Complex d2(0., 0.), e(0., 0.);


	GSL::Vector a(3), b(3), c(3);
	a[0] = 1.; a[1] = 0.; a[2] = 0.;
	b[0] = 0.; b[1] = 1.; b[2] = 0.;
	c[0] = 0.; c[1] = 0.; c[2] = 1.;
	Crystal cr(6*a, 6*b, 6*c);

	std::cout << cr.lat.scale*cr.lat.lat << std::endl;
	size_t Nk = cr.calc_nk(1e-6, sqrt(0.015), lm {4, 0});
	size_t Nr = cr.calc_nr(1e-6, sqrt(0.015), lm {4, 0});
//	std::cout << "Nr = " << Nr << std::endl;
	cr.calc_Rn(Nr/3);
//	std::cout << "Nk = " << Nk << std::endl;
	cr.calc_Kn(Nk/3);


	GSL::Vector tau(3);

	Atom C1(mesh, tau*cr.lat.scale*cr.lat.lat);

	tau[0] = 0.5;
	tau[1] = 0.5;
	tau[2] = 0.;
	kp[0] = 0.;

	Atom C2(mesh, tau*cr.lat.scale*cr.lat.lat);
	tau[0] = 0.;
	tau[1] = 0.5;
	tau[2] = 0.5;
	Atom C3(mesh, tau*cr.lat.scale*cr.lat.lat);
	tau[0] = 0.5;
	tau[1] = 0.;
	tau[2] = 0.5;
	Atom C4(mesh, tau*cr.lat.scale*cr.lat.lat);

	C1.set_Z(6);
	C2.set_Z(6);
	C3.set_Z(6);
	C4.set_Z(6);

	cr.add_atoms(std::vector<Atom> {C1, C2, C3, C4});

	Bloch_sum bs( lm {3, 2}, sqrt(0.015), cr);

	Bloch_summed_structure_constant B(cr, lm {3, 2}, lm {1, -1});

	std::vector<std::vector<Atom>> nn = cr.calc_nearest_neighbours();
	cr.atoms[0].set_MT(nn[0][0].get_pos().norm()/2);
	cr.atoms[1].set_MT(nn[1][0].get_pos().norm()/2);
	cr.atoms[2].set_MT(nn[2][0].get_pos().norm()/2);
	cr.atoms[3].set_MT(nn[3][0].get_pos().norm()/2);

	double at_vol = 0;
	for(size_t idx = 0; idx < cr.atoms.size(); idx++){
		at_vol += GSL::pow_int(cr.atoms[idx].get_MT(), 3);
	}
	at_vol *= 4*M_PI/3;

	for(size_t idx = 0; idx < cr.atoms.size(); idx++){
		cr.atoms[idx].set_AS( std::cbrt(cr.volume/at_vol) *
		cr.atoms[idx].get_MT());
	}
	at_vol = 0.;
	for(size_t idx = 0; idx < cr.atoms.size(); idx++){
		at_vol += GSL::pow_int(cr.atoms[idx].get_AS(), 3);
	}
	at_vol *= 4*M_PI/3;

	std::cout << "Atomic volumes = " << at_vol << " a.u.^3" << std::endl;
	std::cout << "Cell volume = " << cr.volume << " a.u.^3" << std::endl;


	cr.atoms[0].mesh = Logarithmic_mesh(cr.atoms[0].get_AS(), 500);
	cr.atoms[1].mesh = Logarithmic_mesh(cr.atoms[1].get_AS(), 500);
	cr.atoms[2].mesh = Logarithmic_mesh(cr.atoms[2].get_AS(), 500);
	cr.atoms[3].mesh = Logarithmic_mesh(cr.atoms[3].get_AS(), 500);

	double kappa = sqrt(0.015);
	Augmented_spherical_wave aw1s(kappa, 1, lm {0, 0}, UP, cr.atoms[0], cr.atoms);
	Augmented_spherical_wave aw2s(kappa, 2, lm {0, 0}, UP, cr.atoms[0], cr.atoms);
	Augmented_spherical_wave aw2p(kappa, 2, lm {1, 0}, UP, cr.atoms[0], cr.atoms);
	aw1s.core_state = true;
	Potential pot(cr.atoms);
	pot.initial_pot();

	for(size_t idx = 0; idx < cr.atoms.size(); idx++){
		std::cout << "Muffin-tin radius of atom  " << idx << " : "<< cr.atoms[idx].get_MT()
		<< std::endl;
		std::cout << "Atomic sphere radius of atom " << idx << " : " << cr.atoms[1].get_AS()
		<< std::endl;
	}

	aw1s.set_up(pot);
	aw2s.set_up(pot);
	aw2p.set_up(pot);

	GSL::Vector r(3);


	std::ofstream out_file;
	out_file.open("check_ASW.dat");
	out_file << "# r\tV(r)\tr*R1s(r)\tr*R2s(r)\tr*R2p(r)" << std::endl;
	for(size_t i = 0; i < 2000; i++){
		r = 0.001*i*tau*cr.lat.scale*cr.lat.lat;
		out_file << std::setprecision(8) << r.norm() << " " << pot(r)<< " "
		<< aw1s(r) << " " << aw2s(r) << " " << aw2p(r) << std::endl;
	}
	out_file.close();

	return 0;
}
