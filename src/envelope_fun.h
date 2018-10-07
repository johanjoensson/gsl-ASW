#ifndef ENVELOPE_FUN_H
#define ENVELOPE_FUN_H
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"
#include "utils.h"

// Class/Functor for representing envelope funcitons
class Envelope_function{
public:
    GSL::Vector center;
    lm l;
    double kappa;

    Envelope_function();
    Envelope_function(const GSL::Vector center, const lm l, const double kappa);
    double operator()(const GSL::Vector r);
    virtual double barred_fun(double x)
        {return 0.*x;}
};

class Envelope_Hankel : public Envelope_function{
public:
    double barred_fun(double x);
};

class Envelope_Bessel : public Envelope_function{
public:
    double barred_fun(double x);
};

class Envelope_Neumann : public Envelope_function{
public:
    double barred_fun(double x);
};

/*******************************************************************************
* Hankel functions                                                             *
*******************************************************************************/
// One center integral
double atomic_integral(Envelope_Hankel H1, Envelope_Hankel H2);
// Two center integrals
double atomic_integral(Envelope_Hankel H1, Envelope_Bessel J2);
double atomic_integral(Envelope_Bessel J1, Envelope_Hankel H2);

/*******************************************************************************
* Neumann functions                                                            *
*******************************************************************************/
// One center integral
double atomic_integral(Envelope_Neumann N1, Envelope_Neumann N2);
// Two center integrals
double atomic_integral(Envelope_Neumann N1, Envelope_Bessel J2);
double atomic_integral(Envelope_Bessel J1, Envelope_Neumann N2);

// Three center integral
double atomic_integral(Envelope_Bessel J1, Envelope_Bessel J2);

#endif // ENVELOPE_FUN_H
