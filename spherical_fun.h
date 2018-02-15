#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <gsl/gsl_vector.h>

unsigned long int factorial(int n);
double real_spherical_hankel(int l, double x);
double cubic_harmonic(int l, int m, gsl_vector r);

#endif //SPHERICAL_FUN_H
