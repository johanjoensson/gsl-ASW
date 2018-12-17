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
    GSL::Matrix XH1, XS1, XH2, XS2, XH3, XS3;

    void add_states(const Atom& at, const double kappa);
    GSL::Complex H_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
    GSL::Complex S_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
public:
    Simulation();
    Simulation(const Crystal& crystal, const XC_FUN func, const double kappa);

    void set_up_X_matrices();
    void set_up_H(const GSL::Vector& kp);
    void set_up_S(const GSL::Vector& kp);
    void calc_eigen();
};
#endif // SIMULATION_H
