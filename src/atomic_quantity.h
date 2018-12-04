#ifndef ATOMIC_QUANTITY_H
#define ATOMIC_QUANTITY_H
#include <vector>
#include "../../GSL-lib/src/vector.h"
#include "atom.h"
#include "xc_func.h"

class Atomic_quantity{
    friend class Augmented_spherical_wave;

public:
    Atomic_quantity() : sites(), val(){};
    Atomic_quantity(Atomic_quantity&) = default;
    Atomic_quantity(Atomic_quantity&&) = default;
    virtual ~Atomic_quantity() = default;

    Atomic_quantity(const std::vector<Atom>& atoms);

    std::vector<Atom> sites;
    std::vector<std::vector<double>> val;

    double operator()(const GSL::Vector& r);

    Atomic_quantity& operator=(const Atomic_quantity&) = default;
    Atomic_quantity& operator=(Atomic_quantity&&) = default;
};

class Potential : public Atomic_quantity{
    std::vector<std::vector<double>> electrostatic, exchange_correlation;

public:
    Xc_func xc_fun;
    void initial_pot(unsigned int nel, double vol);

    double MT_0 = 0;

    Potential() : Atomic_quantity(), electrostatic(), exchange_correlation(), xc_fun(), MT_0(0){};
    Potential(std::vector<Atom>& atoms);
    ~Potential(){};

    Potential(Potential&) = default;
    Potential(Potential&&) = default;
    Potential& operator=(const Potential&) = default;
    Potential& operator=(Potential&&) = default;

    void set_xc_fun(XC_FUN xc_func);
};

class Density : public Atomic_quantity{
    std::vector<double> valence, core;

public:
    Density() : Atomic_quantity(), valence(), core(){};
    ~Density() {};
    Density(Density&) = default;
    Density(Density&&) = default;

    Density& operator=(const Density&) = default;
    Density& operator=(Density&&) = default;

    Density(std::vector<Atom>& atoms);
};
#endif // ATOMIC_QUANTITY_H
