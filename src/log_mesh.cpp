#include "log_mesh.h"
#include "../../GSL-lib/src/special_functions.h"

Logarithmic_mesh::Logarithmic_mesh()
 : A_p(), B_p(), r(), r2(), drx()
{}

Logarithmic_mesh::Logarithmic_mesh(double radius, size_t num_points)
	: A_p(0.02),
      B_p(radius/(GSL::exp(this->A_p*(static_cast<double>(num_points - 1))).val - 1.)),
      r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	double tmp = 0;

	for(size_t i = 0; i < num_points; i++){
		tmp = B_p*(GSL::exp(A_p*static_cast<double>(i)).val - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A_p*(tmp + this->B_p);
	}
}

Logarithmic_mesh::Logarithmic_mesh(double A_n, double radius,
    size_t num_points)
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
