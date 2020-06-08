#ifndef AUGMENTED_SPHERICAL_WAVE_H
#define AUGMENTED_SPHERICAL_WAVE_H
#include "GSLpp/vector.h"
#include "utils.h"
#include "atomic_quantity.h"
#include "atom.h"
#include "augmented_fun.h"

#include <unordered_set>

class Augmented_spherical_wave{
private:

public:
    Atom center;
    std::vector<Atom> off_centers;
    double kappa;
    size_t index;
    lm l;
    spin s;
    bool core_state = false;

    Augmented_Hankel H;
    std::vector<std::unordered_set<Augmented_Bessel>> J;

    Augmented_spherical_wave(double kappa, size_t index, lm l, spin s,
        const Atom& center, const std::vector<Atom>& off_centers);

    void set_up(Potential &v);

    double operator()(const GSL::Vector &r) const;
};

#endif // AUGMENTED_SPHERICAL_WAVE_H
