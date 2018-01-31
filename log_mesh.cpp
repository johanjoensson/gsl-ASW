#include <iostream>

#include "log_mesh.h"
#include <gsl/gsl_sf_exp.h>

Logarithmic_mesh::Logarithmic_mesh(double radius, uint num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = 0.02;
	this->B = radius/(gsl_sf_exp(this->A*(num_points - 1)) - 1.);
	double tmp;

	for(uint i = 0; i < num_points; i++){
		tmp = B*(gsl_sf_exp(A*i) - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}

Logarithmic_mesh::Logarithmic_mesh(double A, double radius, uint num_points)
	: r(num_points,0), r2(num_points,0), drx(num_points,0)
{
	this->A = A;
	this->B = radius/(gsl_sf_exp(this->A*(num_points - 1)) - 1.);
	double tmp;

	for(uint i = 0; i < num_points; i++){
		tmp = this->B*(gsl_sf_exp(this->A*i) - 1);
		this->r[i] = tmp;
		this->r2[i] = tmp*tmp;
		this->drx[i] = this->A*(tmp + this->B);
	}
}
/*
int main()
{
	uint len = 100;
	Logarithmic_mesh lmesh(0.02,8.0,len);

	for(uint i = 0; i < len; i++)
	{
		std::cout << lmesh.r[i] << "\t"<< lmesh.r2[i] << "\t"
		<< lmesh.drx[i] << std::endl;
	}

	return 0;
}
*/
