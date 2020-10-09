#ifndef ENVELOPE_FUN_H
#define ENVELOPE_FUN_H
#include "GSLpp/vector.h"
#include "GSLpp/special_functions.h"
#include "utils.h"
#include "atom.h"
#include "site.h"

// Class/Functor for representing envelope funcitons
class Envelope_function{
protected:
    lm l_m;
    double kappa_m;

public:
    Envelope_function(lm l_n, double kappa_n)
     : l_m(l_n), kappa_m(kappa_n)
    {}
    Envelope_function(const Envelope_function&) = default;
    Envelope_function(Envelope_function&&) = default;

    Envelope_function& operator=(const Envelope_function&) = default;
    Envelope_function& operator=(Envelope_function&&) = default;
    virtual ~Envelope_function(){};

    double operator()(const GSL::Vector r);
    virtual double barred_fun(const double x) const
        {return 0.*x;}

    double kappa() const {return kappa_m;}
    double& kappa() {return kappa_m;}

    lm l() const {return l_m;}
    lm& l() {return l_m;}

    friend double atomic_integral(const Envelope_function&, const Envelope_function&, const double rs);
};

class Envelope_Hankel;
class Envelope_Bessel;
class Envelope_Neumann;

class Envelope_Hankel : public Envelope_function{
public:
    Envelope_Hankel(lm l_n, double kappa_n)
     : Envelope_function(l_n, kappa_n)
    {};

    double barred_fun(const double x) const override;

    friend double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Hankel& H2, const double rs);
    friend double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);
    friend double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);

};

class Envelope_Bessel : public Envelope_function{
public:
    Envelope_Bessel(lm l_n, double kappa_n)
     : Envelope_function(l_n, kappa_n)
    {};

    double barred_fun(const double x) const override;
    friend double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);
    friend double atomic_integral(const Envelope_Bessel& J1, const Envelope_Bessel& J2, const double rs);
    friend double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);
};

class Envelope_Neumann : public Envelope_function{
public:
    // Envelope_Neumann()
     // :Envelope_function()
     // {};
    Envelope_Neumann(lm l_n, double kappa_n)
     : Envelope_function(l_n, kappa_n)
    {};

    double barred_fun(const double x) const override;

};

/*******************************************************************************
* Hankel functions                                                             *
*******************************************************************************/
// One center integral
double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Hankel& H2, const double rs);
double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);
// Two center integrals
double atomic_integral(const Envelope_Hankel& H1, const Envelope_Bessel& J2, const double rs);
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Hankel& H2, const double rs);

/*******************************************************************************
* Neumann functions                                                            *
*******************************************************************************/
// One center integral
double atomic_integral(const Envelope_Neumann& N1, const Envelope_Neumann& N2, const double rs);
// Two center integrals
double atomic_integral(const Envelope_Neumann& N1, const Envelope_Bessel& J2, const double rs);
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Neumann& N2, const double rs);

// Three center integral
double atomic_integral(const Envelope_Bessel& J1, const Envelope_Bessel& J2, const double rs);

/*******************************************************************************
* Wronskkians of envelope functions (barred spherical functions)               *
*******************************************************************************/
double wronskian(const Envelope_Hankel& f, const Envelope_Hankel& g, const double rs);

#endif // ENVELOPE_FUN_H
