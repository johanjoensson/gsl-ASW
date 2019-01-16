#include "log_mesh.h"
#include "GSLpp/special_functions.h"
#include <stdexcept>

Logarithmic_mesh::Logarithmic_mesh()
 : A_p(), B_p(), r(), r2(), drx()
{}

Logarithmic_mesh::Logarithmic_mesh(double radius,
    size_t num_points, double A_n)
	: A_p(A_n),
      B_p(radius/(GSL::exp(this->A_p*(static_cast<double>(num_points - 1))).val - 1.)),
      r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = this->B_p*(GSL::exp(this->A_p*i).val - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A_p*(tmp + this->B_p);
	}
}

// Trapezoidal rule for integrals on a logarithmic mesh
double Logarithmic_mesh::integrate(std::vector<double>& f)
{
    if(r.size() != f.size()){
    	throw std::runtime_error("Mesh size does not match function size!\n");
    }
    double res = 0;
    for(size_t i = 1; i < r.size(); i++){
    	res += 0.5*(f[i - 1] + f[i])*(r[i] - r[i - 1]);
    }
    return res;
}

// Composite Simpson's rule for integrals on a logarithmic mesh
double Logarithmic_mesh::integrate_simpson(std::vector<double>& f)
{
    if(r.size() != f.size()){
    	throw std::runtime_error("Mesh size does not match function size!\n");
    }
    double A = 0., B = 0., C = 0., res = 0.;//, Ap = 0;
    // Assume even number of mesh-points
    for(size_t i = 1; i < ( r.size() )/2; i++){
    	A = (f[2*i] - f[2*i - 1])/((r[2*i] - r[2*i - 1])*(r[2*i] - r[2*i - 2])) - 
	    (f[2*i - 1] - f[2*i - 2])/((r[2*i - 1] - r[2*i - 2])*(r[2*i] - r[2*i - 2]));
/*	    
	Ap = f[2*i - 2]/((r[2*i - 1] - r[2*i - 2])*(r[2*i] - r[2*i - 2])) - 
	     f[2*i - 1]/((r[2*i - 1] - r[2*i - 2])*(r[2*i] - r[2*i - 1])) + 
	     f[2*i]    /((r[2*i] - r[2*i - 2])*(r[2*i] - r[2*i - 1]));
*/

	B = (f[2*i] - f[2*i - 1])/(r[2*i] - r[2*i - 1]) - A*(r[2*i - 1] + r[2*i]);
	C = f[2*i - 2] - A*r2[2*i - 2] - B*r[2*i - 2];

        res += 	A/3*(r[2*i]*r2[2*i] - r[2*i - 2]*r2[2*i - 2]) + 
		B/2*(r2[2*i] - r2[2*i - 2]) + 
		C*(r[2*i] - r[2*i - 2]);
    }
    return res;
}
