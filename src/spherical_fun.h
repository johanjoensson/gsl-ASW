#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include <gsl/gsl_vector.h>

struct lm {
	int l;
	int m;
};

unsigned long int factorial(int n);
double real_spherical_hankel(lm l, double x);
double cubic_harmonic(lm l, gsl_vector r);


std::ostream& operator << ( std::ostream&, const lm& );  

#endif //SPHERICAL_FUN_H
