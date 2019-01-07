#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "utils.h"
#include "ewald_int.h"
#include "GSLpp/vector.h"
#include "GSLpp/special_functions.h"

// Functor for representing Spherical Bessel functions
class Spherical_function{
public:
    lm l;
    Spherical_function(): l() {};
    Spherical_function(const lm l_n):l(l_n){};
    Spherical_function(const Spherical_function& f) : Spherical_function(f.l)
    {};
    Spherical_function(Spherical_function&& f) : Spherical_function()
    {
        std::swap(l, f.l);
    };

    Spherical_function& operator=(const Spherical_function& f) = default;
    Spherical_function& operator=(Spherical_function&& f) = default;

    virtual ~Spherical_function(){};
    virtual double operator()(const double x)
        {return 0.*x;}
};

double wronskian(Spherical_function& a, Spherical_function& b, double r);

class Hankel_function : public Spherical_function
{
public:
    Hankel_function() : Spherical_function(){};
    Hankel_function(const lm l_n) : Spherical_function(l_n){};
    double operator()(const double x);
};

class Bessel_function : public Spherical_function
{
public:
    Bessel_function() : Spherical_function(){};
    Bessel_function(const lm l_n) : Spherical_function(l_n){};
    double operator()(const double x);
};

class Neumann_function : public Spherical_function
{
public:
    Neumann_function() : Spherical_function(){};
    Neumann_function(const lm l_n) : Spherical_function(l_n){};
    double operator()(const double x);
};

class Integral_Hankel_function : public Hankel_function
{
private:
    Ewald_integral I;
public:
    Integral_Hankel_function() : Hankel_function(), I(){ I.set_kappa(1.0);};
    Integral_Hankel_function(const lm l_n) : Hankel_function(l_n), I(){};
    double operator()(const double x);
    void set_ewald_param(const double eta);
};

unsigned long int factorial(int n);
GSL::Result cubic_harmonic(lm l, const GSL::Vector& r);

std::ostream& operator << ( std::ostream& os, const lm& l);
#endif //SPHERICAL_FUN_H
