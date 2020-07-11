#ifndef ATOMIC_QUANTITY_H
#define ATOMIC_QUANTITY_H
#include <vector>
#include <functional>
#include "GSLpp/vector.h"
#include "atom.h"
#include "crystal.h"
#include "xc_func.h"
#include "site.h"

class Atomic_quantity{
    friend class Augmented_spherical_wave;
protected:

    const Crystal_t<3,Atom>& cr;
    const std::vector<Logarithmic_mesh>& at_meshes;
    std::vector<std::vector<double>> val;

public:
    Atomic_quantity(const Crystal_t<3, Atom>& cr, const std::vector<Logarithmic_mesh>& at_meshes);

    double operator()(const GSL::Vector& r);

    std::vector<double>& sphere(const size_t i){return val[i];}
};

class Potential : public Atomic_quantity{
    std::vector<std::vector<double>> electrostatic, exchange_correlation;
    std::function<double(const size_t, const double)> at_pot;
    double MT_0;

    double Xi0(const Site_t<3>& j, const double r);

public:
    Xc_func xc_fun;
    void initial_pot(double vol);

    double MT0(){return MT_0;};

    Potential(const Crystal_t<3, Atom>& cryst, const std::vector<Logarithmic_mesh>& at_meshes,
        std::function<double(const size_t Z, const double r)> atomic_potential =
        [](const size_t Z, const double r){
            return -2.*static_cast<double>(Z)/r;
        });

    void set_xc_fun(XC_FUN xc_func);
};

class Density : public Atomic_quantity{
    std::vector<double> valence, core;

public:
    Density(const Crystal_t<3,Atom>&, std::vector<Logarithmic_mesh>&);
};
#endif // ATOMIC_QUANTITY_H
