#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "spherical_fun.h"
#include "numerov_solver.h"
#include "structure_const.h"
#include "numerov_solver.h"
#include "atom.h"
#include "crystal.h"

int l;
int z;
double kappa = std::sqrt(0.015);

double v_at(double r)
{
	return -2.*z/r;
}

double v_eff(double r)
{
	return l*(l+1.)/(r*r) + v_at(r);
}



int main()
{
	z = 29;
	gsl_vector *r = gsl_vector_alloc(3*sizeof(double));

	Logarithmic_mesh m(1,1000);
	gsl_vector_set(r, 0, 0.);
	gsl_vector_set(r, 1, 0.5);
	gsl_vector_set(r, 2, 0.);
	Atom a1(1, 1, z, *r, m);
	gsl_vector_set(r, 0, 0.);
	gsl_vector_set(r, 1, 0.2);
	gsl_vector_set(r, 2, 0.);
	Atom a2(1, 1, z, *r, m);

	gsl_vector tmp;
	double x = 0., y = 0., z = 0;
	std::vector<Atom*> atoms {&a1, &a2};
	for(Atom *at : atoms){
		tmp = at->get_pos();
		x = gsl_vector_get(&tmp, 0);
		y = gsl_vector_get(&tmp, 1);
		z = gsl_vector_get(&tmp, 2);
		std::cout << "Atom at position (" << x << ", " << y << ", " << z << ")!" << std::endl;
	}
	gsl_vector_free(r);

	return 0;
}
