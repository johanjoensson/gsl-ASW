#include "log_mesh.h"
#include "GSLpp/special_functions.h"
#include <stdexcept>

Logarithmic_mesh::Logarithmic_mesh(double radius, size_t num_points, double A_n)
	: Mesh(num_points), A_p(A_n),
    B_p(radius/GSL::expm1(A_n*(static_cast<double>(num_points) - 1))),
    drx_p(num_points, 0)
{
	double r;

	for(unsigned int i = 0; i < num_points; i++){
		r = this->B_p*GSL::expm1(this->A_p*i);
		this->x_p[i] = r;
		this->x2_p[i] = GSL::pow_int(r, 2);
		this->drx_p[i] = this->A_p*(r + this->B_p);
	}
}

// Trapezoidal rule for integrals on a logarithmic mesh
double Logarithmic_mesh::integrate(const std::vector<double>& f) const
{
    if(this->size() != f.size()){
    	throw std::runtime_error("Mesh size does not match function size!\n");
    }
    double res = 0;
    for(size_t i = 1; i < this->size(); i++){
    	res += 0.5*(f[i - 1] + f[i])*drx_p[i-1];
    }
    return res;
}

// Composite Simpson's rule for integrals on a logarithmic mesh
double Logarithmic_mesh::integrate_simpson(const std::vector<double>& f) const
{
    if(this->size() != f.size()){
    	throw std::runtime_error("Mesh size does not match function size!\n");
    }
    double res = 0.;
    for(size_t i = 1; i < ( this->size())/2; i++){
        res += f[2*i - 2]*this->drx_p[2*i - 2] +
               4*f[2*i - 1]*this->drx_p[2*i - 1] +
			   f[2*i]*this->drx_p[2*i];
    }
	if(this->size() % 2 != 0){
		size_t i = this->size();
		res += f[i - 3]*this->drx_p[i - 3] + 4*f[i - 2]*this->drx_p[i - 2] +
		       f[i - 1]*this->drx_p[i - 1];
	}
    return 1./3 * res;
}
