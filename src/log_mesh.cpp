#include "log_mesh.h"
#include "spherical_fun.h"
#include "../../GSL-lib/src/special_functions.h"

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

Logarithmic_mesh::Logarithmic_mesh(double A, double radius, unsigned int num_points)
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
