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
		Mesh() : x_p(), x2_p(), dx_p() {}
		Mesh(const Mesh&) = default;
		Mesh(Mesh&&) = default;

		virtual ~Mesh() = default;

		Mesh& operator=(const Mesh&) = default;
		Mesh& operator=(Mesh&&) = default;
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
		std::vector<double>::const_iterator x_begin() const {return x_p.begin();}
		std::vector<double>::iterator x_end() {return x_p.end();}
		std::vector<double>::const_iterator x_end() const {return x_p.end();}
		std::vector<double>::reverse_iterator x_rbegin() {return x_p.rbegin();};
		std::vector<double>::const_reverse_iterator x_rbegin() const {return x_p.rbegin();};
		std::vector<double>::reverse_iterator x_rend() {return x_p.rend();};
		std::vector<double>::const_reverse_iterator x_rend() const {return x_p.rend();};
		double x_back() const {return x_p.back();}
		double x2(const size_t i) const {return x2_p[i];}
		std::vector<double>::iterator x2_begin() {return x2_p.begin();}
		std::vector<double>::const_iterator x2_begin() const {return x2_p.begin();}
		std::vector<double>::iterator x2_end() {return x2_p.end();}
		std::vector<double>::const_iterator x2_end() const {return x2_p.end();}
		std::vector<double>::reverse_iterator x2_rbegin() {return x2_p.rbegin();};
		std::vector<double>::const_reverse_iterator x2_rbegin() const {return x2_p.rbegin();};
		std::vector<double>::reverse_iterator x2_rend() {return x2_p.rend();};
		std::vector<double>::const_reverse_iterator x2_rend() const {return x2_p.rend();};
		double x2_back() const {return x2_p.back();}
		double dx() const {return dx_p;}

		size_t size() const {return x_p.size();}

		virtual double integrate(const std::vector<double>& f) const = 0;
		virtual double integrate_simpson(const std::vector<double>& f) const = 0;
};

class Logarithmic_mesh : public Mesh {
		double A_p;
		double B_p;
		std::vector<double> drx_p;
	public:
		Logarithmic_mesh() : Mesh(), A_p(), B_p(), drx_p() {}
		Logarithmic_mesh(const Logarithmic_mesh&) = default;
		Logarithmic_mesh(Logarithmic_mesh&&) = default;

		~Logarithmic_mesh() override = default;

		Logarithmic_mesh(double radius, size_t num_points, double A = 0.01);

		Logarithmic_mesh& operator=(const Logarithmic_mesh&) = default;
		Logarithmic_mesh& operator=(Logarithmic_mesh&&) = default;

		double integrate(const std::vector<double>& f) const override;
		double integrate_simpson(const std::vector<double>& f) const override;

		double A() const {return A_p;};
		double B() const {return B_p;};
		double r(const size_t i) const {return x(i);};
		std::vector<double>::iterator r_begin() {return x_begin();};
		std::vector<double>::const_iterator r_begin() const {return x_begin();};
		std::vector<double>::iterator r_end() {return x_end();};
		std::vector<double>::const_iterator r_end() const {return x_end();};
		std::vector<double>::reverse_iterator r_rbegin() {return x_rbegin();};
		std::vector<double>::const_reverse_iterator r_rbegin() const {return x_rbegin();};
		double r_back() const {return x_back();};
		double r2(const size_t i) const {return x2(i);};
		std::vector<double>::iterator r2_begin() {return x2_begin();};
		std::vector<double>::const_iterator r2_begin() const {return x2_begin();};
		std::vector<double>::iterator r2_end() {return x2_end();};
		std::vector<double>::const_iterator r2_end() const {return x2_end();};
		double r2_back() const {return x2_back();};
		double drx(const size_t i) const {return drx_p[i];};
		std::vector<double>::iterator drx_begin() {return drx_p.begin();};
		std::vector<double>::const_iterator drx_begin() const {return drx_p.begin();};
		std::vector<double>::iterator drx_end() {return drx_p.end();};
		std::vector<double>::const_iterator drx_end() const {return drx_p.end();};
		std::vector<double>::reverse_iterator drx_rbegin() {return drx_p.rbegin();};
		std::vector<double>::const_reverse_iterator drx_rbegin() const {return drx_p.rbegin();};
		std::vector<double>::reverse_iterator drx_rend() {return drx_p.rend();};
		std::vector<double>::const_reverse_iterator drx_rend() const {return drx_p.rend();};
		double drx_back() const {return drx_p.back();};
};
#endif //LOG_MESH_H
