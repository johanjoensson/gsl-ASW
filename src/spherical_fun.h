#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include "../../GSL-lib/src/vector.h"

struct lm {
	int l;
	int m;
};

unsigned long int factorial(int n);
double real_spherical_hankel(lm l, double x);
double cubic_harmonic(lm l, GSL::Vector& r);

double exp_gsl(double x);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
