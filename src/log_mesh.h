#ifndef LOG_MESH_H
#define LOG_MESH_H

#include <cstddef>
#include <vector>

class Logarithmic_mesh {
		double A_p;
		double B_p;
		std::vector<double> r_p;
		std::vector<double> r2_p;
		std::vector<double> drx_p;
	public:
		Logarithmic_mesh();
		Logarithmic_mesh(double radius, size_t num_points, double A = 0.02);

		double integrate(std::vector<double>& f);
		double integrate_simpson(std::vector<double>& f);

		double A() const {return A_p;};
		double B() const {return B_p;};
		double r(const size_t i) const {return r_p[i];};
		double r_back() const {return r_p.back();};
		double r2(const size_t i) const {return r2_p[i];};
		double r2_back() const {return r2_p.back();};
		double drx(const size_t i) const {return drx_p[i];};
		double drx_back() const {return drx_p.back();};
		size_t size() const {return r_p.size();};
};
#endif //LOG_MESH_H
