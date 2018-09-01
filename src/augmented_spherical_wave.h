#ifndef AUGMENTED_SPHERICAL_WAVE_H
#define AUGMENTED_SPHERICAL_WAVE_H
#include "../../GSL-lib/src/vector.h"
#include "utils.h"
#include "atomic_quantity.h"
#include "atom.h"
#include "augmented_fun.h"

#include <unordered_set>

/***************************************************************************//**
* Class for representing an augmented spherical wave.\n
* Contains:\n
* __center__ - Atomic site the wave is centered on, on-center augmentation sphere.\n
* __off_centers__ - Atomic sites where the wave is augmented in off-center spheres.\n
* __kappa_ - Energy of waves in interstitial region.\n
* __n__ - principal quantum number of the wave.\n
* __l__ - Combined orbital and magnetic quantum number.\n
* __s__ - Spin direction.\n
* __core_state_ - Is this wave describing a core state or not? Core states have
* no augmentation in off-center spheres.\n
* __H__ - On center, augmented Hankel function.\n
* __J__ - Off-center, augmented Bessel functions (unsorted).\n
*******************************************************************************/
class Augmented_spherical_wave{
private:

public:
    Atom center;
    std::vector<Atom> off_centers;
    double kappa;
    unsigned int n;
    lm l;
    spin s;
    bool core_state = false;

    Augmented_Hankel H;
    std::vector<std::unordered_set<Augmented_Bessel>> J;

    Augmented_spherical_wave();
    Augmented_spherical_wave(double kappa, unsigned int n, lm l, spin s,
        Atom center, std::vector<Atom> off_centers);

    //! Initialise wave using potential v.
    void set_up(Potential &v);

    //! Obtain value of wave at point r.
    double operator()(const GSL::Vector &r);
};

#endif // AUGMENTED_SPHERICAL_WAVE_H
