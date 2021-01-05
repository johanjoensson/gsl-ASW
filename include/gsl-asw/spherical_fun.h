#ifndef SPHERICAL_FUN_H
#define SPHERICAL_FUN_H
#include <ostream>
#include "utils.h"
#include "ewald_int.h"
#include "GSLpp/vector.h"
#include "GSLpp/special_functions.h"

// Functor for representing Spherical Bessel functions
class Spherical_function{
protected:
    lm l_m;
public:
    Spherical_function() : l_m() {};
    Spherical_function(const lm l_n):l_m(l_n){}
    Spherical_function(const Spherical_function& f) = default;
    Spherical_function(Spherical_function&& f) = default;

    Spherical_function& operator=(const Spherical_function& f) = default;
    Spherical_function& operator=(Spherical_function&& f) = default;

    virtual ~Spherical_function(){}
    virtual double operator()(const double x) const
        {return 0.*x;}

    void set_l(const lm l_n)
        {l_m = l_n;};
    lm l() const
        {return l_m;};
};

double wronskian(Spherical_function a, Spherical_function b, double r);

class Hankel_function : public Spherical_function
{
public:
    Hankel_function() : Spherical_function(){}
    Hankel_function(const lm l_n) : Spherical_function(l_n){}
    double operator()(const double x) const override;
};

class Bessel_function : public Spherical_function
{
public:
    Bessel_function() : Spherical_function(){};
    Bessel_function(const lm l_n) : Spherical_function(l_n){};
    double operator()(const double x) const override;
};

class Neumann_function : public Spherical_function
{
public:
    Neumann_function() : Spherical_function(){};
    Neumann_function(const lm l_n) : Spherical_function(l_n){};
    double operator()(const double x) const override;
};

class Integral_Hankel_function : public Hankel_function
{
private:
    Ewald_integral I;
public:
    Integral_Hankel_function() : Hankel_function(), I(){};
    using Hankel_function::operator();
    double operator()(const double kappa, const double eta, const double x) const;
};

unsigned long int factorial(int n);

GSL::Result cubic_harmonic(const lm& l, const GSL::Vector& r);
GSL::Result cubic_harmonic(const int l, const int m, const double cos_theta, const double phi);

#endif //SPHERICAL_FUN_H
