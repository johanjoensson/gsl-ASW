#ifndef ENVELOPE_FUN_H
#define ENVELOPE_FUN_H
#include "../../GSL-lib/src/vector.h"
#include "../../GSL-lib/src/special_functions.h"
#include "utils.h"
#include "atom.h"

// Class/Functor for representing envelope funcitons
class Envelope_function{
public:
    Atom center;
    lm l;
    double kappa;

    Envelope_function();
    Envelope_function(const Atom& center, lm l, double kappa);

    /*
    Envelope_function(Envelope_function &a);
    Envelope_function(Envelope_function &&a);

    Envelope_function& operator=(const Envelope_function& a);
    Envelope_function& operator=(Envelope_function&& a);
    */
    double operator()(const GSL::Vector r);
    virtual double barred_fun(const double x)
        {return 0.*x;}
};

class Envelope_Hankel : public Envelope_function{
public:
    Envelope_Hankel()
     :Envelope_function()
     {};
    Envelope_Hankel(const Atom& center, lm l, double kappa)
     : Envelope_function(center, l, kappa)
     {};

     /*
    Envelope_Hankel(Envelope_Hankel& a)
     : Envelope_function(a)
     {}
    Envelope_Hankel(Envelope_Hankel&& a)
     : Envelope_function(a)
     {}
     */
    double barred_fun(const double x) const;
};

class Envelope_Bessel : public Envelope_function{
public:
    Envelope_Bessel()
     :Envelope_function()
     {};
    Envelope_Bessel(const Atom& center, lm l, double kappa)
     : Envelope_function(center, l, kappa)
     {};
    /*
    Envelope_Bessel(Envelope_Bessel& a)
     : Envelope_function(a)
     {}
    Envelope_Bessel(Envelope_Bessel&& a)
     : Envelope_function(a)
     {}
    */
    double barred_fun(const double x) const;
};

class Envelope_Neumann : public Envelope_function{
public:
    Envelope_Neumann()
     :Envelope_function()
     {};
    Envelope_Neumann(const Atom& center, lm l, double kappa)
     : Envelope_function(center, l, kappa)
     {};
    /*
    Envelope_Neumann(Envelope_Neumann& a)
     : Envelope_function(a)
     {}
    Envelope_Neumann(Envelope_Neumann&& a)
     : Envelope_function(a)
     {}
    */
    double barred_fun(const double x) const;
};

/*******************************************************************************
* Hankel functions                                                             *
*******************************************************************************/
// One center integral
double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Hankel& H2);
// Two center integrals
double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2);
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Hankel& H2);

/*******************************************************************************
* Neumann functions                                                            *
*******************************************************************************/
// One center integral
double atomic_integral(const Envelope_Neumann& N1, const Envelope_Neumann& N2);
// Two center integrals
double atomic_integral(const Envelope_Neumann& N1, const Envelope_Bessel& J2);
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Neumann& N2);

// Three center integral
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Bessel& J2);

#endif // ENVELOPE_FUN_H
