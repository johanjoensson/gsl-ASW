#ifndef NUMEROV_SOLVER_H
#define NUMEROV_SOLVER_H

#include <vector>
#include "log_mesh.h"

class Numerov_solver{
		double (*v_ext) (double r);
	public:
		Numerov_solver();

		std::vector<double> v;
		void set_v_ext(double (*new_v_ext) (double r));
		std::vector<double> solve_left(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E);
		std::vector<double> solve_right(Logarithmic_mesh &mesh, int i_lr, std::vector<double> &init_cond, double E);

};

#endif //NUMEROV_SOLVER_H
