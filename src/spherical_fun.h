#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

struct lm {
	int l;
	int m;
};

unsigned long int factorial(int n);
double real_spherical_hankel(lm l, double x);
double cubic_harmonic(lm l, gsl_vector r);

double exp_gsl(double x);


std::ostream& operator << ( std::ostream&, const lm& );
std::ostream& operator << ( std::ostream&, const gsl_complex& );
gsl_complex operator + (gsl_complex a, gsl_complex b);
gsl_complex operator - (gsl_complex a, gsl_complex b);
gsl_complex operator * (gsl_complex a, gsl_complex b);
gsl_complex operator / (gsl_complex a, gsl_complex b);

gsl_complex operator + (gsl_complex a, double b);
gsl_complex operator - (gsl_complex a, double b);
gsl_complex operator * (gsl_complex a, double b);
gsl_complex operator / (gsl_complex a, double b);

gsl_complex operator + (double a, gsl_complex b);
gsl_complex operator - (double a, gsl_complex b);
gsl_complex operator * (double a, gsl_complex b);
gsl_complex operator / (double a, gsl_complex b);

gsl_vector& operator + (gsl_vector& a, gsl_vector& b);
#endif //SPHERICAL_FUN_H
