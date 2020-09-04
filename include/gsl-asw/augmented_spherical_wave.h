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
    // std::vector<Hankel_container>& Hs_m;
    // std::vector<Bessel_container>& Bs_m;
    Site_t<3> center_m;
    double kappa_m;
    lm l_m;
    spin s_m;
    bool core_state_m = false;

public:
    Augmented_spherical_wave() = default;
    Augmented_spherical_wave(const Augmented_spherical_wave&) = default;
    Augmented_spherical_wave(Augmented_spherical_wave&&) = default;

    // Augmented_spherical_wave(std::vector<Hankel_container>& Hs_n, std::vector<Bessel_container>& Bs_n, Site_t<3> center, double kappa, lm l, spin s);
    Augmented_spherical_wave(Site_t<3> center, double kappa, lm l, spin s, bool core_state = false)
     : center_m(center), kappa_m(kappa), l_m(l), s_m(s), core_state_m(core_state)
    {}

    ~Augmented_spherical_wave() = default;

    Augmented_spherical_wave& operator=(const Augmented_spherical_wave&) = default;
    Augmented_spherical_wave& operator=(Augmented_spherical_wave&&) = default;

    Site_t<3>& center(){return center_m;}
    Site_t<3> center() const {return center_m;}
    double kappa() const {return kappa_m;}
    lm l() const {return l_m;}
    spin s() const {return s_m;}
    bool core_state() const {return core_state_m;}

    // double operator()(const GSL::Vector &r) const;
};

#endif // AUGMENTED_SPHERICAL_WAVE_H
