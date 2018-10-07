#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "utils.h"
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"

// Functor for representing Spherical Bessel functions
class Spherical_function{
public:
    lm l;
    Spherical_function();
    Spherical_function(const lm l);
    virtual double operator()(const double x)
        {return 0.*x;}
};

double wronskian(Spherical_function a, Spherical_function b, double r);

class Hankel_function : public Spherical_function
{
public:
    Hankel_function() : Spherical_function(){};
    Hankel_function(const lm l) : Spherical_function(l){};
    double operator()(const double x);
};

class Bessel_function : public Spherical_function
{
public:
    Bessel_function() : Spherical_function(){};
    Bessel_function(const lm l) : Spherical_function(l){};
    double operator()(const double x);
};

class Neumann_function : public Spherical_function
{
public:
    Neumann_function() : Spherical_function(){};
    Neumann_function(const lm l) : Spherical_function(l){};
    double operator()(const double x);
};

unsigned long int factorial(int n);
GSL::Result cubic_harmonic(lm l, const GSL::Vector& r);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
