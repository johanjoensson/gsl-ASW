#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <thread>
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
#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"
#include "GSLpp/complex.h"
#include "GSLpp/special_functions.h"
#include "GSLpp/error.h"
#include "augmented_spherical_wave.h"
#include "atomic_quantity.h"
#include "envelope_fun.h"
#include "xc_func.h"
#include "k-mesh.h"

void k_iteration(const GSL::Vector& kp, Simulation& sim)
{
	sim.calc_eigen(kp);
}

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


	double kappa = std::sqrt(0.015);


	GSL::Vector a = {0.0, 0.5, 0.5}, b = {0.5, 0.0, 0.5}, c = {0.5, 0.5, 0.0};
	// GSL::Vector a = {1.0, 0.0, 0.0}, b = {0.0, 1.0, 0.0}, c = {0.0, 0.0, 1.0};
	std::cout << a << "\n";
	std::cout << b << "\n";
	std::cout << c << "\n";
	// Crystal_t<3, Atom> cr(Lattice_t<3>({16*a, 16*b, 16*c}));
	Crystal_t<3, Atom> cr(Lattice_t<3>({3*a, 3*b, 3*c}));
	for(auto row : cr.lat().recip_lat()){
		std::cout << row <<"\n";
	}
	std::cout << "\n";

	std::cout << "Crystal volume = " << cr.volume() << " (a.u.)^3\n";

	std::cout << "Calculating Rmax" << "\n";
	double Rmax = calc_Rmax(cr.volume(), kappa, lm {5, 0}, 5e-14);
	std::cout << "Calculating Kmax" << "\n";
	double Kmax = calc_Kmax(cr.volume(), kappa, lm {5, 0}, 5e-14);
	std::cout << "Rmax = " << Rmax << ", Kmax = " << Kmax << "\n";


	std::cout << "Setting R-vectors" << std::endl;
	cr.set_Rn(Rmax);
	std::cout << "Setting K-vectors" << std::endl;
	cr.set_Kn(Kmax);
	std::cout << "Number of R-vectors " << cr.Rn_vecs().size() << "\n";
	std::cout << "Number of K-vectors " << cr.Kn_vecs().size() << "\n";

	std::cout << "Setting up atoms" << std::endl;
	GSL::Vector tau(3);
	Atom C1(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());

	tau[0] = 0.25;
	tau[1] = 0.25;
	tau[2] = 0.25;
	Atom C2(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());

	tau[0] = 0.;
	tau[1] = 0.5;
	tau[2] = 0.5;
	Atom C3(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());

	tau[0] = 0.5;
	tau[1] = 0.;
	tau[2] = 0.5;
	Atom C4(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());
	Atom C5(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());
	Atom C6(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());
	Atom C7(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());
	Atom C8(/*Logarithmic_mesh(),*/ tau*cr.lat().lat());

	C1.set_Z(21);
	C2.set_Z(21);
	C3.set_Z(6);
	C4.set_Z(6);
	C5.set_Z(6);
	C6.set_Z(6);
	C7.set_Z(6);
	C8.set_Z(6);

	cr.set_size({1, 1, 1});
	cr.add_basis({C1});//, C2});//, C3, C4, C5, C6, C7, C8});
	// cr.add_sites({{0, 0, 0}, {0.5, 0, 0}, {0, 0.5, 0}, {0, 0, 0.5}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}, {0.5, 0.5, 0.5}});
	cr.add_sites({{0, 0, 0}, {0.25, 0.25, 0.25}});

	std::cout << "Crystal contains " << cr.sites().size() << " sites\n";
	std::cout << "Crystal contains " << cr.atoms().size() << " inequivalent atoms\n";

	Simulation sim(cr, LDA, {kappa});

	std::cout << "Setting up X matrices\n" << std::flush;
	sim.set_up_X_matrices();
	std::cout << "Initialising K-mesh" << std::endl;
	K_mesh kmesh(cr.lat().recip_lat());


/*
	std::cout << "Testing Gaunt coefficients\n";
	for(lm l1 = {2, 1, -1}; l1 != lm {3, 0, 0}; l1++){
		for(lm l2 = {2, 1, -1}; l2 != lm {3, 0, 0}; l2++){
			for(lm l3 = {4, 0, 0}; l3 != lm {5, 0, 0}; l3++){
				std::cout << "l1 = " << l1 << "\n";
				std::cout << "l2 = " << l2 << "\n";
				std::cout << "l3 = " << l3 << "\n";
				std::cout << "\tGaunt coefficient = " << gaunt(l1, l2, l3) << "\n";
			}
		}
	}
*/

	// kmesh.generate_mesh(10, 10, 10);
	// kmesh.generate_mesh({{0,0,0}, {1,1,1}}, 1);
	kmesh.generate_mesh({{0.5, 0.5, 0.5}, {0, 0, 0}, {0, 0.5, 0.5}, {0, 0, 0}, {0, 0.5, 0.5}, {0, 0, 0}, {0.5, 0, 0}}, 100);
	std::cout << "K-mesh generated" << std::endl;

/*
	size_t n_threads = 1;
	std::vector<std::thread> thread_pool(n_threads);
*/

    std::ofstream of;
	of.open("canonical_bands.dat", std::ios::trunc);
    for(const auto& kp :  kmesh.k_points){
		of << kp[0] << " " << kp[1] << " " << kp[2] << " ";
		for(lm l = {3, 2, -2}; l != lm {4, 0, 0}; l++){
			of << sim.canonical_band(l, kappa, UP, kp) << " ";
		}
		of << "\n";
	}
	of.close();

	// return 0;

	std::cout << "Loop ove k-points\n" << std::endl;
	for(const auto& kp : kmesh.k_points){
/*
		size_t t_id = 0;
		while(true){
			if(!thread_pool[t_id].joinable()){
				thread_pool[t_id] = std::thread(
					[&sim] (const GSL::Vector& k)
					{
*/
						sim.calc_eigen(kp);
/*
					}, kp);
					break;
			}else{
				thread_pool[t_id].join();
			}
			t_id = (t_id + 1) % n_threads;
		}
*/
	}
	sim.print_eigvals(kmesh);
/*
	for(auto&& thread : thread_pool){
		if(thread.joinable()){
			thread.join();
		}
	}

	std::cout << "All threads joined succesfully!\n" << std::endl;
*/
	return 0;
}
