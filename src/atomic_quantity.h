#ifndef ATOMIC_QUANTITY_H
#define ATOMIC_QUANTITY_H
#include <vector>
#include "../../GSL-lib/src/vector.h"
#include "atom.h"

class Atomic_quantity{
    friend class Augmented_spherical_wave;

public:
    Atomic_quantity();
    Atomic_quantity(Atomic_quantity&) = default;
    Atomic_quantity(Atomic_quantity&&) = default;
    virtual ~Atomic_quantity() = default;

    Atomic_quantity(const std::vector<Atom> atoms);

    std::vector<Atom> sites;
    std::vector<std::vector<double>> val;

    double operator()(const GSL::Vector& r);

    Atomic_quantity& operator=(const Atomic_quantity&) = default;
    Atomic_quantity& operator=(Atomic_quantity&&) = default;
};

class Potential : public Atomic_quantity{
    std::vector<std::vector<double>> electrostatic, exchange_correlation;

public:
    void initial_pot();

    double MT_0;

    Potential(std::vector<Atom>& atoms);
    ~Potential(){};

    Potential(Potential&) = default;
    Potential(Potential&&) = default;
    Potential& operator=(const Potential&) = default;
    Potential& operator=(Potential&&) = default;

};

class Density : public Atomic_quantity{
    std::vector<double> valence, core;

public:
    ~Density() {};
    Density(Density&) = default;
    Density(Density&&) = default;

    Density& operator=(const Density&) = default;
    Density& operator=(Density&&) = default;

    Density(std::vector<Atom>& atoms);
};
#endif // ATOMIC_QUANTITY_H
