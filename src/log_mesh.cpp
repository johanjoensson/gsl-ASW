#include "log_mesh.h"
#include "GSLpp/special_functions.h"
#include <stdexcept>

Logarithmic_mesh::Logarithmic_mesh()
 : A_p(), B_p(), r_p(), r2_p(), drx_p()
{}

Logarithmic_mesh::Logarithmic_mesh(double radius, size_t num_points, double A_n)
	: A_p(A_n), B_p(radius/(GSL::exp(A_n*(static_cast<double>(num_points) - 1)).val - 1.)),
      r_p(num_points, 0), r2_p(num_points, 0), drx_p(num_points, 0)
{
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = this->B_p*(GSL::exp(this->A_p*i).val - 1);
		this->r_p[i] = tmp;
		this->r2_p[i] = tmp*tmp;
		this->drx_p[i] = this->A_p*(tmp + this->B_p);
		// this->drx_p[i] = (tmp + this->B_p)*(GSL::exp(A_p).val - 1);
	}
    // this->drx_p[this->size() - 1] = this->A_p*(this->r_back() + this->B_p);
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
    	// res += 0.5*(f[i - 1] + f[i])*(r_p[i] - r_p[i - 1]);
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
