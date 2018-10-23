#include "log_mesh.h"
#include "../../GSL-lib/src/special_functions.h"

Logarithmic_mesh::Logarithmic_mesh()
 : B(), r(), r2(), drx(), A()
{}

Logarithmic_mesh::Logarithmic_mesh(double radius, unsigned int num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = 0.02;
	this->B = radius/(GSL::exp(this->A*(num_points - 1)).val - 1.);
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = B*(GSL::exp(A*i).val - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}

Logarithmic_mesh::Logarithmic_mesh(double A, double radius,
    unsigned int num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = A;
	this->B = radius/(GSL::exp(this->A*(num_points - 1)).val - 1.);
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = this->B*(GSL::exp(this->A*i).val - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}

// Composite Simpson's rule for radial integrals on a logarithmic mesh
double Logarithmic_mesh::radial_integral(std::vector<double>& f)
{
    double t0 = 0., t1 = 0., t2 = 0., res = 0.;
    for(size_t i = 1; i < (r.size() - 1)/2; i++){
        // i is an integer, 2*i has to be even.
        // This part is 2*i even
        t2 = f[2*i - 1]*r2[2*i - 1]*drx[2*i - 2];
        t1 = f[2*i]*r2[2*i - 1]*drx[2*i - 1];
        t0 = f[2*i + 1]*r2[2*i + 1]*drx[2*i];
        res += t2 + 4*t1 + t0;
        // We need to make sure we include the cases where 2*i _should
        // be odd.
        // This part is 2*i odd
        t2 = f[2*i]*r2[2*i]*drx[2*i - 1];
        t1 = f[2*i]*r2[2*i]*drx[2*i];
        t0 = f[2*i + 2]*r2[2*i + 2]*drx[2*i + 1];
        res += t2 + 4*t1 + t0;
    }
    return 1./3 * res;
}
