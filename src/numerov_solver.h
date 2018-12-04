#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <vector>
#include "log_mesh.h"

class Numerov_solver{
		unsigned int find_inversion_point(Logarithmic_mesh &mesh,
			std::vector<double> &v, double e_trial);
		double variational_energy_correction(Logarithmic_mesh &mesh,
			std::vector<double> &v, std::vector<double> &fun, int i_inv,
			double e_trial);
	public:
		Numerov_solver();

		std::vector<double> solve_left(Logarithmic_mesh &mesh,
			std::vector<double> &v, int i_lr, std::vector<double> &init_cond,
			double E);
		std::vector<double> solve_right(Logarithmic_mesh &mesh,
			std::vector<double> &v, int i_lr, std::vector<double> &init_cond,
			double E);

		std::vector<double> solve(Logarithmic_mesh &mesh,
			std::vector<double> &v, std::vector<double> &l_init,
			std::vector<double> &r_init, double& en, int n_nodes);
};

#endif //NUMEROV_SOLVER_H
