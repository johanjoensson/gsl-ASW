#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "utils.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"

unsigned long int factorial(int n);
GSL::Result real_spherical_hankel(lm l, double x);
GSL::Result cubic_harmonic(lm l, const GSL::Vector& r);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
