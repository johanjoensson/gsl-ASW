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

		Logarithmic_mesh(double radius, uint num_points);
		Logarithmic_mesh(double A, double radius, uint num_points);
};
#endif //LOG_MESH_H
