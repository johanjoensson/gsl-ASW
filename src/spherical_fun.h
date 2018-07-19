#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"

struct lm {
	int l;
	int m;
};

unsigned long int factorial(int n);
GSL::Result real_spherical_hankel(lm l, double x);
GSL::Result cubic_harmonic(lm l, GSL::Vector& r);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
