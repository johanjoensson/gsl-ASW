#include "log_mesh.h"
#include "spherical_fun.h"

Logarithmic_mesh::Logarithmic_mesh(double radius, unsigned int num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = 0.02;
	this->B = radius/(exp_gsl(this->A*(num_points - 1)) - 1.);
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = B*(exp_gsl(A*i) - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}

Logarithmic_mesh::Logarithmic_mesh(double A, double radius, unsigned int num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = A;
	this->B = radius/(exp_gsl(this->A*(num_points - 1)) - 1.);
	double tmp;

	for(unsigned int i = 0; i < num_points; i++){
		tmp = this->B*(exp_gsl(this->A*i) - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}
