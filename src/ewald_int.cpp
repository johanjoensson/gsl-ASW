#include "ewald_int.h"
#include <cmath>

Integral_hankel::Integral_hankel(Logarithmic_mesh &mesh)
    : mesh(mesh)
{
    this->kappa = sqrt(0.015);
    this->ewald_param = 1;
}

void Integral_hankel::set_ewald_param(double eta)
{
	this->ewald_param = eta;
}

void Integral_hankel::set_kappa(double kappa)
{
	this->kappa = kappa;
}
