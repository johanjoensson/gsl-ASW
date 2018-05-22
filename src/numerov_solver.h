#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <vector>
#include "log_mesh.h"

/***************************************************************************//**
* A class for solving the radial Schrödinger equation using Numerov's method\n
* Find approximate energy values using variational method ensuring the solution
* has the correct number of nodes\n
*******************************************************************************/
* *class Numerov_solver{
		//! Effective potential
		double (*v_eff) (double r);
		//! Potential inside atomic sphere
		double (*v_at) (double r);
		//! Function for locating point where left and right integrations should meet
		int find_inversion_point(Logarithmic_mesh &mesh, double e_trial);
		//! Approximate correction to energy parameter using variational method
		double variational_energy_correction(Logarithmic_mesh &mesh, std::vector<double> &fun, int i_inv, double e_trial);
	public:
		Numerov_solver();

		//! Set effective potential
		void set_v_eff(double (*new_v_eff) (double r));
		//! Set atomic potential
		void set_v_at(double (*new_v_at) (double r));
		//! Solve radial Schrödinger equation starting from the outer edge of the sphere ending at the inflection point
		std::vector<double> solve_left(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E);
		//! Solve radial Schrödinger equation starting from the leftmost point ending at the inflection point
		std::vector<double> solve_right(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E);

		//! Combine solutions from the left and the right and make sure they
		//! have the correct number of nodes and are continuous and smooth at
		//! the inflection point
		std::vector<double> solve(Logarithmic_mesh &mesh, std::vector<double> &l_init, std::vector<double> &r_init, double &en, int n_nodes);
};

#endif //NUMEROV_SOLVER_H
