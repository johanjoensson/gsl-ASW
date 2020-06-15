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
    std::vector<Hankel_container>& Hs_m;
    std::vector<Bessel_container>& Bs_m;

public:
    Site_t<3> center;
    double kappa;
    lm l;
    spin s;
    bool core_state = false;

    Augmented_spherical_wave(std::vector<Hankel_container>& Hs_n, std::vector<Bessel_container>& Bs_n, Site_t<3> center, double kappa, lm l, spin s);

    double operator()(const GSL::Vector &r) const;
};

#endif // AUGMENTED_SPHERICAL_WAVE_H
