#ifndef LOG_MESH_H
#define LOG_MESH_H

#include <vector>

/***************************************************************************//**
* A class for representing logarithmic meshes\n
* Contains:\n
* __B__ - Parameter determining shape of mesh, calculated from __num_points__ and __A__\n
* __r__ - r-values contained in mesh\n
* __r2__ - r^2-values contained in mesh\n
* __drx__ - Derivative dr/dx evaluated in mesh\n
* __A__ - Parameter controlling spacing between points in mesh\n
*******************************************************************************/
class Logarithmic_mesh {
		double B;
	public:
		std::vector<double> r;
		std::vector<double> r2;
		std::vector<double> drx;

		double A;

		Logarithmic_mesh(double radius, unsigned int num_points);
		Logarithmic_mesh(double A, double radius, unsigned int num_points);
};
#endif //LOG_MESH_H
