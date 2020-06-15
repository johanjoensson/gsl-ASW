#ifndef ENVELOPE_FUN_H
#define ENVELOPE_FUN_H
#include "GSLpp/vector.h"
#include "GSLpp/special_functions.h"
#include "utils.h"
#include "atom.h"

// Class/Functor for representing envelope funcitons
class Envelope_function{
public:
    size_t center;
    lm l;
    double kappa;

    // Envelope_function() : center(), l(), kappa() {};
    Envelope_function(const size_t& center, lm l, double kappa);
    Envelope_function(const Envelope_function&) = default;
    Envelope_function(Envelope_function&&) = default;

    Envelope_function& operator=(const Envelope_function&) = default;
    Envelope_function& operator=(Envelope_function&&) = default;
    virtual ~Envelope_function(){};

    double operator()(const GSL::Vector r);
    virtual double barred_fun(const double x) const
        {return 0.*x;}
};

class Envelope_Hankel : public Envelope_function{
public:
    // Envelope_Hankel()
    //  :Envelope_function()
    //  {};
    Envelope_Hankel(const size_t& center_n, lm l_n, double kappa_n)
     : Envelope_function(center_n, l_n, kappa_n)
     {};

    double barred_fun(const double x) const override;
};

class Envelope_Bessel : public Envelope_function{
public:
    // Envelope_Bessel()
    //  :Envelope_function()
    //  {};
    Envelope_Bessel(const size_t& center_n, lm l_n, double kappa_n)
     : Envelope_function(center_n, l_n, kappa_n)
     {};

    double barred_fun(const double x) const override;
};

class Envelope_Neumann : public Envelope_function{
public:
    // Envelope_Neumann()
     // :Envelope_function()
     // {};
    Envelope_Neumann(const size_t& center_n, lm l_n, double kappa_n)
     : Envelope_function(center_n, l_n, kappa_n)
     {};

    double barred_fun(const double x) const override;
};

/*******************************************************************************
* Hankel functions                                                             *
*******************************************************************************/
// One center integral
double off_atomic_integral(const Envelope_Hankel& H1, const Envelope_Hankel& H2, const double rs);
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

#endif // ENVELOPE_FUN_H
