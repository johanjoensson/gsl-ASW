#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <map>
#include <functional>
#include "atomic_quantity.h"
#include "augmented_spherical_wave.h"
#include "crystal.h"
#include "xc_func.h"
#include "k-mesh.h"

struct GSLVecCompare{
    bool operator()(const GSL::Vector& lhs, const GSL::Vector& rhs) const
    {
        return (lhs.norm<double>() < rhs.norm<double>());
    }
};

class Simulation{
    Crystal cryst;
    Potential pot;
    Density n;
    std::vector<Augmented_spherical_wave> basis_valence;
    std::vector<Augmented_spherical_wave> basis_core;
    std::map<GSL::Vector, GSL::Vector, GSLVecCompare> k_eigenvals;
    GSL::Matrix_cx H, S;
    GSL::Matrix XH1, XS1, XH2, XS2, XH3, XS3;

    void add_states(const Atom& at, const double kappa);
    GSL::Complex H_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
    GSL::Complex S_element(const size_t i1, const size_t i2, const GSL::Vector& kp);
public:
    Simulation();
    Simulation(const Crystal& crystal, const XC_FUN func, const double kappa,
    std::function<double(const size_t, const double)> at_pot =
        [](const size_t Z, const double r)
        {
            return -2.*static_cast<double>(Z)/r;
        });

    void set_up_X_matrices();
    void set_up_H(const GSL::Vector& kp);
    void set_up_S(const GSL::Vector& kp);
    std::pair<GSL::Matrix_cx, GSL::Vector> calc_eigen();

    void add_eigvals(const GSL::Vector& kp, const GSL::Vector& eigvals);
    void print_eigvals(K_mesh& k_mesh);
};
#endif // SIMULATION_H
