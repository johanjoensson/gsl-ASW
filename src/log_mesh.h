#ifndef LOG_MESH_H
#define LOG_MESH_H

#include <vector>

class Logarithmic_mesh {
		double B;
	public:
		std::vector<double> r;
		std::vector<double> r2;
		std::vector<double> drx;

		double A;

		Logarithmic_mesh();
		Logarithmic_mesh(double radius, unsigned int num_points);
		Logarithmic_mesh(double A, double radius, unsigned int num_points);

		double radial_integral(std::vector<double> f);
};
#endif //LOG_MESH_H
