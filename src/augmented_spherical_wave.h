#ifndef AUGMENTED_SPHERICAL_WAVE_H
#define AUGMENTED_SPHERICAL_WAVE_H
#include "../../GSL-lib/src/vector.h"
#include "spherical_fun.h"

enum spin{
    UP,
    DOWN,
    NON_COLLINEAR
};

class Augmented_spherical_wave{
private:
    GSL::Vector on_center;
    double kappa, r_MT, r_AS;
    lm l;
    spin s;

public:
    Augmented_spherical_wave();
};

#endif // AUGMENTED_SPHERICAL_WAVE_H
