#include "log_mesh.h"
#include "GSLpp/special_functions.h"
#include <stdexcept>

Logarithmic_mesh::Logarithmic_mesh(double radius, size_t num_points, double A_n)
	: Mesh(num_points), A_p(A_n),
    B_p(radius/(GSL::exp(A_n*(static_cast<double>(num_points) - 1)).val - 1.)),
    drx_p(num_points, 0)
{
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = this->B_p*(GSL::exp(this->A_p*i).val - 1);
		this->x_p[i] = tmp;
		this->x2_p[i] = tmp*tmp;
		this->drx_p[i] = this->A_p*(tmp + this->B_p);
	}
}

// Trapezoidal rule for integrals on a logarithmic mesh
double Logarithmic_mesh::integrate(std::vector<double>& f)
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
double Logarithmic_mesh::integrate_simpson(std::vector<double>& f)
{
    if(this->size() != f.size()){
    	throw std::runtime_error("Mesh size does not match function size!\n");
    }
    double res = 0.;
    for(size_t i = 1; i < ( this->size())/2; i++){
        res += f[2*i - 2]*this->drx_p[2*i - 2] +
               4*f[2*i - 1]*this->drx_p[2*i - 1] + f[2*i]*this->drx_p[2*i];
    }
    return 1./3 * res;
}
