#ifndef LOG_MESH_H
#define LOG_MESH_H

#include <cstddef>
#include <vector>

class Logarithmic_mesh {
		double A_p;
		double B_p;
	public:
		std::vector<double> r;
		std::vector<double> r2;
		std::vector<double> drx;


		Logarithmic_mesh();
		Logarithmic_mesh(double radius, size_t num_points, double A = 0.02);

		double integrate(std::vector<double>& f);
		double integrate_simpson(std::vector<double>& f);

		double A() const {return A_p;};
		double B() const {return B_p;};
};
#endif //LOG_MESH_H
