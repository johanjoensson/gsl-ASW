#ifndef ATOMIC_QUANTITY_H
#define ATOMIC_QUANTITY_H
#include <vector>
#include <functional>
#include "GSLpp/vector.h"
#include "atom.h"
#include "crystal.h"
#include "xc_func.h"

class Atomic_quantity{
    friend class Augmented_spherical_wave;

public:
    Atomic_quantity() : sites(), val(){}

    Atomic_quantity(const Crystal_t<3, Atom>& cr);

    std::vector<Atom> sites;
    std::vector<std::vector<double>> val;

    double operator()(const GSL::Vector& r);
};

class Potential : public Atomic_quantity{
    std::vector<std::vector<double>> electrostatic, exchange_correlation;
    std::function<double(const size_t, const double)> at_pot;
    double MT_0;

public:
    Xc_func xc_fun;
    void initial_pot(size_t nel, double vol);

    double MT0(){return MT_0;};

    Potential() : Atomic_quantity(), electrostatic(), exchange_correlation(),
        at_pot([](const size_t Z, const double r){
            return 0.*(static_cast<double>(Z) + r);}), MT_0(),
        xc_fun(){}
    Potential(const Crystal_t<3, Atom>& cryst,
        std::function<double(const size_t Z, const double r)> atomic_potential =
        [](const size_t Z, const double r){
            return -2.*static_cast<double>(Z)/r;
        });

    void set_xc_fun(XC_FUN xc_func);
};

class Density : public Atomic_quantity{
    std::vector<double> valence, core;

public:
    Density() : Atomic_quantity(), valence(), core(){}

    Density(const Crystal_t<3,Atom>&);
};
#endif // ATOMIC_QUANTITY_H
