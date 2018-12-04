#ifndef SIMULATION_H
#define SIMULATION_H
#include "../../GSL-lib/src/complex_matrix.h"
#include "atomic_quantity.h"
#include "augmented_spherical_wave.h"
#include "crystal.h"
#include "xc_func.h"

class Simulation{
    Crystal cryst;
    Potential pot;
    Density n;
    std::vector<Augmented_spherical_wave> basis_valence;
    std::vector<Augmented_spherical_wave> basis_core;
    GSL::Complex_Matrix H, S;

    void add_states(Atom& at, double kappa);
    GSL::Complex H_element(const Augmented_spherical_wave& w1, const Augmented_spherical_wave& w2, const GSL::Vector& kp);
    GSL::Complex S_element(const Augmented_spherical_wave& w1, const Augmented_spherical_wave& w2, const GSL::Vector& kp);
public:
    Simulation();
    Simulation(Crystal& crystal, XC_FUN func, double kappa);

    void set_up_H(const GSL::Vector& kp);
    void set_up_S(const GSL::Vector& kp);
    void calc_eigen();
};
#endif // SIMULATION_H
