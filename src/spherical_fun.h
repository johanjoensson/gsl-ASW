#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "utils.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"

//! Somewhat useful implementation of factorial function.
unsigned long int factorial(int n);
//! Implementation of spherical hankel function for imaginary argument x.
GSL::Result real_spherical_hankel(lm l, double x);
//! Implementation of real spherical harmonics (cubic harmonics).
GSL::Result cubic_harmonic(lm l, const GSL::Vector& r);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
