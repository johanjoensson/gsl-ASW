#ifndef LOG_MESH_H
#define LOG_MESH_H

#include <GSLpp/special_functions.h>
#include <cstddef>
#include <vector>

class Mesh{
protected:
		std::vector<double> x_p;
		std::vector<double> x2_p;
		double dx_p;
	public:
		virtual ~Mesh(){}
		Mesh(double x0, double x1, double dx)
		 : x_p(), x2_p(), dx_p(dx)
		{
			for(double x = x0; x <= x1; x += dx){
				x_p.push_back(x);
				x2_p.push_back(GSL::pow_int(x, 2));
			}
		}

		Mesh(size_t n) : x_p(n, 0), x2_p(n, 0), dx_p(0){};
		Mesh(double x0, double x1, size_t n) : Mesh(x0, x1, (x1 - x0)/static_cast<double>(n)) {}

		double x(const size_t i) const {return x_p[i];}
		std::vector<double>::iterator x_begin() {return x_p.begin();}
		std::vector<double>::reverse_iterator x_rbegin() {return x_p.rbegin();};
		double x_back() const {return x_p.back();}
		double x2(const size_t i) const {return x2_p[i];}
		std::vector<double>::iterator x2_begin() {return x2_p.begin();}
		std::vector<double>::reverse_iterator x2_rbegin() {return x2_p.rbegin();};
		double x2_back() const {return x2_p.back();}
		double dx() const {return dx_p;}

		size_t size() const {return x_p.size();}

		double integrate(std::vector<double>& f);
		double integrate_simpson(std::vector<double>& f);
};

class Logarithmic_mesh : public Mesh {
		double A_p;
		double B_p;
		std::vector<double> drx_p;
	public:
		~Logarithmic_mesh() = default;
		Logarithmic_mesh(double radius, size_t num_points, double A = 0.01);

		double integrate(std::vector<double>& f);
		double integrate_simpson(std::vector<double>& f);

		double A() const {return A_p;};
		double B() const {return B_p;};
		double r(const size_t i) const {return x(i);};
		std::vector<double>::iterator r_begin() {return x_begin();};
		std::vector<double>::reverse_iterator r_rbegin() {return x_rbegin();};
		double r_back() const {return x_back();};
		double r2(const size_t i) const {return x2(i);};
		std::vector<double>::iterator r2_begin() {return x2_begin();};
		std::vector<double>::reverse_iterator r2_rbegin() {return r2_rbegin();};
		double r2_back() const {return x2_back();};
		double drx(const size_t i) const {return drx_p[i];};
		std::vector<double>::iterator drx_begin() {return drx_p.begin();};
		std::vector<double>::reverse_iterator drx_rbegin() {return drx_p.rbegin();};
		double drx_back() const {return drx_p.back();};
};
#endif //LOG_MESH_H
